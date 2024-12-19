#[derive(Debug)]
pub struct LweCiphertext<M> {
    pub(crate) a: Vec<M>,
    pub(crate) b: M,
}
