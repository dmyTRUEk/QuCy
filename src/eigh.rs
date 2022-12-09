/// Function for solution of generalized Eigenvalues and Eigenvectors problem.


use nalgebra::{Complex, SymmetricEigen};

use crate::{
    math::{float, Matrix},
    swap,
};



pub trait ExtensionMatrixIsSymetric {
    fn is_symetric(&self) -> bool;
}
impl ExtensionMatrixIsSymetric  for Matrix {
    fn is_symetric(&self) -> bool {
        // self == &self.transpose()
        const THRESHOLD: float = 1e-7;
        for i in 0..self.shape().0 {
            for j in 0..self.shape().1 {
                if i == j { continue; }
                // if self[(i,j)] != self[(j,i)] {
                if ((self[(i,j)] - self[(j,i)])/self[(i,j)]).abs() > THRESHOLD {
                    return false;
                }
            }
        }
        true
    }
}

/// Translated from: https://github.com/scipy/scipy/blob/a00a93c3fee03cd26ebf7b4cba58c2a11d8f52c2/scipy/linalg/_decomp.py#L270-L620
#[allow(non_snake_case)]
pub fn eigh(A: &Matrix, B: &Matrix, _verbose: bool) -> (Vec<Complex<float>>, Matrix) {
    let verbose = false;
    if verbose {
        println!("A:\n{}", A);
        println!("B:\n{}", B);
    }

    // FOR TEST: override it
    // let A = Matrix2::<f64>::new(1.0, 2.0, 2.0, 1.0);
    // let B = Matrix2::<f64>::new(3.0, 1.0, 1.0, 3.0);

    // TODO: enable this!
    assert!(A.is_symetric(), "{A}");
    assert!(B.is_symetric());

    let ed = SymmetricEigen::new(B.clone()); // TODO: optimize: remove `.clone()`
    let (Phi_B, Lambda_B) = (ed.eigenvectors, ed.eigenvalues);
    // let Phi_B = Phi_B.transpose();
    if verbose {
        println!("Phi_B:\n{}", Phi_B);
        println!("Lambda_B:\n{}", Lambda_B);
    }
    // let lambda_B = Lambda_B.diagonal()
    let a = Lambda_B.map(|x| x.sqrt());
    if verbose {
        println!("a:\n{}", a);
    }
    // let a = np.nan_to_num(a) + 0.0001;
    // Lambda_B_squareRoot = np.diag(lambda_B**0.5);

    // let Lambda_B_squareRoot = Matrix2::from_diagonal(&a_);
    // let Phi_B_hat = Phi_B * Lambda_B_squareRoot.try_inverse().unwrap(); // TODO: check
    let Lambda_B_squareRoot: Matrix = Matrix::from_diagonal(&a);
    if verbose {
        println!("Lambda_B_squareRoot:\n{}", Lambda_B_squareRoot);
        println!("Lambda_B_squareRoot.inv:\n{}", Lambda_B_squareRoot.clone().try_inverse().unwrap());
    }
    let Phi_B_hat: Matrix = Phi_B * Lambda_B_squareRoot.try_inverse().unwrap(); // TODO: check
    if verbose {
        println!("Phi_B_hat:\n{}", &Phi_B_hat);
    }

    let A_hat = Phi_B_hat.transpose() * A * &Phi_B_hat;
    if verbose {
        println!("A_hat:\n{}", A_hat);
    }
    let ed = SymmetricEigen::new(A_hat);
    let (Phi_A, Lambda_A) = (ed.eigenvectors, ed.eigenvalues);
    if verbose {
        println!("Phi_A:\n{}", Phi_A);
        println!("Lambda_A:\n{}", Lambda_A);
    }
    let Lambda = Lambda_A;
    let Phi = Phi_B_hat * Phi_A;
    if verbose {
        println!("Phi:\n{Phi}");
    }
    let Lambda: Vec<Complex<float>> = Lambda.iter().map(|&x| Complex::new(x, 0.0)).collect();
    if verbose {
        println!("Lambda:\n{:?}", Lambda);
        println!("\n");
    }
    // todo!("SORT EIG VALS (Lambda) & EIG VECS (Phi) by descending order");

    // assert!(Lambda.shape());
    // assert!(Phi.shape());
    assert!(Lambda.len() == Phi.shape().0);
    assert!(Lambda.len() == Phi.shape().1);

    // sort
    let mut Lambda = Lambda;
    let mut Phi = Phi;
    for i in 0..Lambda.len() {
        for j in i+1..Lambda.len() {
            if Lambda[i].re > Lambda[j].re {
                swap!(Lambda[i], Lambda[j]);
                for k in 0..Lambda.len() {
                    swap!(Phi[(k,i)], Phi[(k,j)]);
                }
            }
        }
    }
    // dbg!(&Lambda);

    return (Lambda, Phi);
}

