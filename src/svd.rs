use nalgebra::{Const, DVector, Dyn, Matrix, SMatrix, SVector, VecStorage};

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

pub fn jacobi_decomposition<const R: usize>(
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

fn calculate_left_singv<const S: usize, const R: usize>(
    a_mat: &SMatrix<f64, S, R>,
    sing_values: &DVector<f64>,
    v_mat: &Matrix<f64, Const<R>, Dyn, VecStorage<f64, Const<R>, Dyn>>,
) -> Matrix<f64, Const<S>, Dyn, VecStorage<f64, Const<S>, Dyn>> {
    let mut column_vec = Vec::with_capacity(sing_values.nrows());
    for i in 0..sing_values.nrows() {
        column_vec.push((a_mat * v_mat.column(i)) / sing_values[i]);
    }
    let u_mat = Matrix::<f64, Const<S>, Dyn, _>::from_columns(&column_vec);
    return u_mat;
}

pub fn svd_decomp<const S: usize, const R: usize>(
    a_mat: &SMatrix<f64, S, R>,
    min_value: f64,
) -> (
    Matrix<f64, Const<S>, Dyn, VecStorage<f64, Const<S>, Dyn>>,
    DVector<f64>,
    Matrix<f64, Const<R>, Dyn, VecStorage<f64, Const<R>, Dyn>>,
) {
    let at_a = a_mat.transpose() * a_mat;

    let (eig_vals, eig_vecs) = jacobi_decomposition::<R>(&at_a, 1e-12, 10000);
    let mut sing_pairs = eig_vals
        .data
        .as_slice()
        .iter()
        .map(|a| a.sqrt())
        .zip(eig_vecs.column_iter().map(|a| a.clone_owned()))
        .filter(|pair| pair.0 > min_value)
        .collect::<Vec<(f64, SMatrix<f64, R, 1>)>>();

    sing_pairs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

    let sing_values =
        DVector::from_iterator(sing_pairs.len(), sing_pairs.iter().map(|pair| pair.0));

    let compact_v = Matrix::<f64, Const<R>, Dyn, _>::from_columns(
        &sing_pairs
            .iter()
            .map(|pair| pair.1)
            .collect::<Vec<Matrix<f64, Const<R>, Const<1>, _>>>(),
    );

    let compact_u = calculate_left_singv(&a_mat, &sing_values, &compact_v);

    (compact_u, sing_values, compact_v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Matrix2x1, Matrix2x3, Matrix3, UniformNorm, Vector3};

    #[test]
    fn simple_max() {
        let mat_a = Matrix3::new(94., 42., 65., 42., 54., 23., 65., 23., 79.);
        let mat_b = Matrix3::new(57., 93., 26., 93., 78., 15., 26., 15., 62.);

        assert_eq!(max_lower_element::<3>(&mat_a), (65., (2, 0)));
        assert_eq!(max_lower_element::<3>(&mat_b), (93., (1, 0)));
    }

    #[test]
    fn test_jacobi() {
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

    #[test]
    fn test_svd() {
        let mat_a = Matrix2x3::new(12.2, 84.3, 48.5, 64.7, 59.2, 46.4);
        let (u, s, v) = svd_decomp(&mat_a, 0.0001);
        assert_eq!(u.shape(), (2, 2));
        assert_eq!(s.shape(), (2, 1));
        assert_eq!(v.shape(), (3, 2));

        let sing_values = Matrix2x1::new(133.255, 41.1639);
        assert!((sing_values - s).apply_norm(&UniformNorm) < 1e-3);
    }
}
