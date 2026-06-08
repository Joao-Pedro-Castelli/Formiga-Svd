pub struct Matrix {
    mat_a: Box<[f64]>,
    size: (usize, usize), // number of lines and number of columns
}

impl Matrix {
    pub fn new(n: usize, m: usize) -> Self {
        let mut mat_a: Vec<f64> = Vec::with_capacity(n * m);
        mat_a.resize(n * m, 0.);

        Matrix {
            mat_a: mat_a.into_boxed_slice(),
            size: (n, m),
        }
    }
}
