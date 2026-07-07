use nalgebra::Matrix2x3;

mod svd;

fn main() {
    let a_mat = Matrix2x3::new(12.2, 84.3, 48.5, 64.7, 59.2, 46.4);
    let (u, s, v) = svd::svd_decomp(&a_mat, 0.05);
    println!("U = {}", u);
    println!("Sigma = {}", s);
    println!("V = {}", v);
    return ();
}
