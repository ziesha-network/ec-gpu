use std::fmt::Write;
use std::mem;

use ec_gpu::GpuField;

static COMMON_SRC: &str = include_str!("cl/common.cl");
static FIELD_SRC: &str = include_str!("cl/field.cl");
static FIELD2_SRC: &str = include_str!("cl/field2.cl");
static EC_SRC: &str = include_str!("cl/ec.cl");

/// Generates the source for FFT and Multiexp operations.
pub fn gen_base_code<Fr: GpuField, Fq: GpuField, L: Limb>() -> String {
    vec![common(), gen_ec_source::<Fr, Fq, L>()].join("\n\n")
}

/// Returns CUDA/OpenCL source-code that contains definitions/functions that are shared across
/// fields.
///
/// It needs to be called before any other function like [`field`] or [`gen_ec_source`] is called,
/// as it contains deinitions, used in those.
pub fn common() -> String {
    COMMON_SRC.to_string()
}

/// Generates the source for the elliptic curve and group operations, as defined by `E`.
///
/// The code from the [`common()`] call needs to be included before this on is used.
pub fn gen_ec_source<Fr: GpuField, Fq: GpuField, L: Limb>() -> String {
    vec![
        field::<Fr, L>("Fr"),
        field::<Fq, L>("Fq"),
        field2("Fq2", "Fq"),
        ec("Fq", "G1"),
        ec("Fq2", "G2"),
    ]
    .join("\n\n")
}

fn field2(field2: &str, field: &str) -> String {
    String::from(FIELD2_SRC)
        .replace("FIELD2", field2)
        .replace("FIELD", field)
}

fn ec(field: &str, point: &str) -> String {
    String::from(EC_SRC)
        .replace("FIELD", field)
        .replace("POINT", point)
}

pub fn field<F, L: Limb>(name: &str) -> String
where
    F: GpuField,
{
    [
        params::<F, L>(),
        field_add_sub_nvidia::<F, L>().expect("preallocated"),
        String::from(FIELD_SRC),
    ]
    .join("\n")
    .replace("FIELD", name)
}

fn const_field<L: Limb>(name: &str, limbs: Vec<L>) -> String {
    format!(
        "CONSTANT FIELD {} = {{ {{ {} }} }};",
        name,
        limbs
            .iter()
            .map(|l| l.value().to_string())
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn params<F, L: Limb>() -> String
where
    F: GpuField,
{
    let one = L::one_limbs::<F>(); // Get Montgomery form of F::one()
    let p = L::modulus_limbs::<F>(); // Get field modulus in non-Montgomery form
    let r2 = L::calculate_r2::<F>();
    let limbs = one.len(); // Number of limbs
    let inv = L::calc_inv(p[0]);
    let limb_def = format!("#define FIELD_limb {}", L::opencl_type());
    let limbs_def = format!("#define FIELD_LIMBS {}", limbs);
    let limb_bits_def = format!("#define FIELD_LIMB_BITS {}", L::bits());
    let p_def = const_field("FIELD_P", p);
    let r2_def = const_field("FIELD_R2", r2);
    let one_def = const_field("FIELD_ONE", one);
    let zero_def = const_field("FIELD_ZERO", vec![L::zero(); limbs]);
    let inv_def = format!("#define FIELD_INV {}", inv.value());
    let typedef = "typedef struct { FIELD_limb val[FIELD_LIMBS]; } FIELD;".to_string();
    [
        limb_def,
        limbs_def,
        limb_bits_def,
        inv_def,
        typedef,
        one_def,
        p_def,
        r2_def,
        zero_def,
    ]
    .join("\n")
}

/// Trait to implement limbs of different underlying bit sizes.
pub trait Limb: Sized + Clone + Copy {
    /// The underlying size of the limb, e.g. `u32`
    type LimbType: Clone + std::fmt::Display;
    /// Returns the value representing zero.
    fn zero() -> Self;
    /// Returns a new limb.
    fn new(val: Self::LimbType) -> Self;
    /// Returns the raw value of the limb.
    fn value(&self) -> Self::LimbType;
    /// Returns the bit size of the limb.
    fn bits() -> usize {
        mem::size_of::<Self::LimbType>() * 8
    }
    /// Returns a tuple with the strings that PTX is using to describe the type and the register.
    fn ptx_info() -> (&'static str, &'static str);
    /// Returns the type that OpenCL is using to represent the limb.
    fn opencl_type() -> &'static str;
    /// Returns the limbs that represent the multiplicative identity of the given field.
    fn one_limbs<F: GpuField>() -> Vec<Self>;
    /// Returns the field modulus in non-Montgomery form as a vector of `Self::LimbType` (least
    /// significant limb first).
    fn modulus_limbs<F: GpuField>() -> Vec<Self>;
    /// Calculate the `INV` parameter of Montgomery reduction algorithm for 32/64bit limbs
    /// * `a` - Is the first limb of modulus.
    fn calc_inv(a: Self) -> Self;
    /// Returns the limbs that represent `R ^ 2 mod P`.
    fn calculate_r2<F: GpuField>() -> Vec<Self>;
}

/// A 32-bit limb.
#[derive(Clone, Copy)]
pub struct Limb32(u32);
impl Limb for Limb32 {
    type LimbType = u32;
    fn zero() -> Self {
        Self(0)
    }
    fn new(val: Self::LimbType) -> Self {
        Self(val)
    }
    fn value(&self) -> Self::LimbType {
        self.0
    }
    fn ptx_info() -> (&'static str, &'static str) {
        ("u32", "r")
    }
    fn opencl_type() -> &'static str {
        "uint"
    }
    fn one_limbs<F: GpuField>() -> Vec<Self> {
        F::one().into_iter().map(Self::new).collect()
    }
    fn modulus_limbs<F: GpuField>() -> Vec<Self> {
        F::modulus().into_iter().map(Self::new).collect()
    }
    fn calc_inv(a: Self) -> Self {
        let mut inv = 1u32;
        for _ in 0..31 {
            inv = inv.wrapping_mul(inv);
            inv = inv.wrapping_mul(a.value());
        }
        Self(inv.wrapping_neg())
    }
    fn calculate_r2<F: GpuField>() -> Vec<Self> {
        F::r2().into_iter().map(Self::new).collect()
    }
}

/// A 64-bit limb.
#[derive(Clone, Copy)]
pub struct Limb64(u64);
impl Limb for Limb64 {
    type LimbType = u64;
    fn zero() -> Self {
        Self(0)
    }
    fn new(val: Self::LimbType) -> Self {
        Self(val)
    }
    fn value(&self) -> Self::LimbType {
        self.0
    }
    fn ptx_info() -> (&'static str, &'static str) {
        ("u64", "l")
    }
    fn opencl_type() -> &'static str {
        "ulong"
    }
    fn one_limbs<F: GpuField>() -> Vec<Self> {
        F::one()
            .chunks(2)
            .map(|chunk| Self::new(((chunk[1] as u64) << 32) + (chunk[0] as u64)))
            .collect()
    }

    fn modulus_limbs<F: GpuField>() -> Vec<Self> {
        F::modulus()
            .chunks(2)
            .map(|chunk| Self::new(((chunk[1] as u64) << 32) + (chunk[0] as u64)))
            .collect()
    }

    fn calc_inv(a: Self) -> Self {
        let mut inv = 1u64;
        for _ in 0..63 {
            inv = inv.wrapping_mul(inv);
            inv = inv.wrapping_mul(a.value());
        }
        Self(inv.wrapping_neg())
    }
    fn calculate_r2<F: GpuField>() -> Vec<Self> {
        F::r2()
            .chunks(2)
            .map(|chunk| Self::new(((chunk[1] as u64) << 32) + (chunk[0] as u64)))
            .collect()
    }
}

/// Generates PTX-Assembly implementation of FIELD_add_/FIELD_sub_
fn field_add_sub_nvidia<F, L: Limb>() -> Result<String, std::fmt::Error>
where
    F: GpuField,
{
    let mut result = String::new();
    let (ptx_type, ptx_reg) = L::ptx_info();

    writeln!(result, "#if defined(OPENCL_NVIDIA) || defined(CUDA)\n")?;
    for op in &["sub", "add"] {
        let len = L::one_limbs::<F>().len();

        writeln!(
            result,
            "DEVICE FIELD FIELD_{}_nvidia(FIELD a, FIELD b) {{",
            op
        )?;
        if len > 1 {
            write!(result, "asm(")?;
            writeln!(result, "\"{}.cc.{} %0, %0, %{};\\r\\n\"", op, ptx_type, len)?;

            for i in 1..len - 1 {
                writeln!(
                    result,
                    "\"{}c.cc.{} %{}, %{}, %{};\\r\\n\"",
                    op,
                    ptx_type,
                    i,
                    i,
                    len + i
                )?;
            }
            writeln!(
                result,
                "\"{}c.{} %{}, %{}, %{};\\r\\n\"",
                op,
                ptx_type,
                len - 1,
                len - 1,
                2 * len - 1
            )?;

            write!(result, ":")?;
            for n in 0..len {
                write!(result, "\"+{}\"(a.val[{}])", ptx_reg, n)?;
                if n != len - 1 {
                    write!(result, ", ")?;
                }
            }

            write!(result, "\n:")?;
            for n in 0..len {
                write!(result, "\"{}\"(b.val[{}])", ptx_reg, n)?;
                if n != len - 1 {
                    write!(result, ", ")?;
                }
            }
            writeln!(result, ");")?;
        }
        writeln!(result, "return a;\n}}")?;
    }
    writeln!(result, "#endif")?;

    Ok(result)
}
