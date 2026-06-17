use nalgebra::SMatrix;

fn find_max_symmetric<const R: usize>(a_mat: &SMatrix<f64, R, R>) -> (f64, (usize, usize)) {
    let mut max = a_mat[(0, 1)].abs();
    let mut index = (0, 1);

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

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Matrix3;

    #[test]
    fn simple_max() {
        let mat_a = Matrix3::new(94., 42., 65., 42., 54., 23., 65., 23., 79.);

        assert_eq!(find_max_symmetric::<3>(&mat_a), (65., (2, 0)));
    }
}
