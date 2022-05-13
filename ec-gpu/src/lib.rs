/// Describes how to generate the gpu sources for a Field.
pub trait GpuField {
    /// Returns `1` as a vector of 32bit limbs.
    fn one() -> Vec<u32>;

    /// Returns `R ^ 2 mod P` as a vector of 32bit limbs.
    fn r2() -> Vec<u32>;

    /// Returns the field modulus in non-Montgomery form (least significant limb first).
    fn modulus() -> Vec<u32>;
}
