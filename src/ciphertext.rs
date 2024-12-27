use crate::poly::FFTPoly;

#[derive(Debug)]
pub struct LweCiphertext<M> {
    pub(crate) a: Vec<M>,
    pub(crate) b: M,
}

#[derive(Debug)]
pub struct NtruScalarCiphertext {
    pub(crate) ct: FFTPoly,
}

// TODO: approximate decomposition
#[derive(Debug)]
pub struct NtruVectorCiphertext {
    pub(crate) vec: Vec<FFTPoly>,
}
