use nalgebra::{SMatrix, SVector};

fn max_lower_element<const R: usize>(a_mat: &SMatrix<f64, R, R>) -> (f64, (usize, usize)) {
    let mut max = a_mat[(1, 0)].abs();
    let mut index = (1, 0);

    for (i, row) in a_mat.view((1, 0), (R - 1, R - 1)).row_iter().enumerate() {
        for (j, element) in row.columns(0, i + 1).column_iter().enumerate() {
            if element[(0, 0)].abs() > max {
                max = element[(0, 0)].abs();
                index = (i + 1, j);
            }
        }
    }

    (max, index)
}

fn calculate_trigonometric<const R: usize>(
    a_mat: &SMatrix<f64, R, R>,
    p: usize,
    q: usize,
) -> (f64, f64) {
    let big_phi = (a_mat[(q, q)] - a_mat[(p, p)]) / (2. * a_mat[(p, q)]);
    let tang = match big_phi.abs() > 1e-15 {
        true => 1. / (big_phi + big_phi.signum() * (big_phi.powi(2) + 1.).sqrt()),
        false => 1.,
    };

    let cosi = 1. / (tang.powi(2) + 1.).sqrt();
    let sino = tang * cosi;

    (cosi, sino)
}

fn jacobi_decomposition<const R: usize>(
    a_mat: &SMatrix<f64, R, R>,
    tol: f64,
    kmax: usize,
) -> (SVector<f64, R>, SMatrix<f64, R, R>) {
    let mut mat_ak = a_mat.clone();
    let mut mat_v = SMatrix::<f64, R, R>::identity();

    for _ in 0..kmax {
        let (max, index) = max_lower_element(&mat_ak);
        if max < tol {
            break;
        }

        let (c, s) = calculate_trigonometric(&mat_ak, index.0, index.1);

        let mut mat_u = SMatrix::identity();
        mat_u[(index.0, index.0)] = c;
        mat_u[(index.0, index.1)] = s;
        mat_u[(index.1, index.0)] = -s;
        mat_u[(index.1, index.1)] = c;

        let tmp = mat_v * mat_u;
        mat_v = tmp;

        mat_ak = mat_v.transpose() * a_mat * mat_v;
    }

    (mat_ak.diagonal(), mat_v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Matrix3, UniformNorm, Vector3};

    #[test]
    fn simple_max() {
        let mat_a = Matrix3::new(94., 42., 65., 42., 54., 23., 65., 23., 79.);
        let mat_b = Matrix3::new(57., 93., 26., 93., 78., 15., 26., 15., 62.);

        assert_eq!(max_lower_element::<3>(&mat_a), (65., (2, 0)));
        assert_eq!(max_lower_element::<3>(&mat_b), (93., (1, 0)));
    }

    #[test]
    fn test_svd() {
        let mat_a = Matrix3::new(94., 42., 65., 42., 54., 23., 65., 23., 79.);
        let mat_b = Matrix3::new(57., 93., 26., 93., 78., 15., 26., 15., 62.);

        // TODO: test eigenvectors
        // let eig_vec_a = Matrix3::new();

        // vector needs to be sorted
        let eig_value_a = Vector3::new(15.833, 40.4584, 170.7086);
        let (mut calc_eigv_a, _) = jacobi_decomposition::<3>(&mat_a, 1e-12, 10000);
        calc_eigv_a
            .as_mut_slice()
            .sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((eig_value_a - calc_eigv_a).apply_norm(&UniformNorm) < 1e-3);

        // let eig_vec_b = Matrix3::new();

        let eig_value_b = Vector3::new(-27.1315, 55.3938, 168.7377);
        let (mut calc_eigv_b, _) = jacobi_decomposition::<3>(&mat_b, 0.00001, 1000);
        calc_eigv_b
            .as_mut_slice()
            .sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((eig_value_b - calc_eigv_b).apply_norm(&UniformNorm) < 1e-3);
    }
}
