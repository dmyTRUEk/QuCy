use std::ops::{Index, IndexMut};

/// Main calculations

use float_extras::f64::erf;
use nalgebra::{DMatrix, Complex};

use crate::{
    math::{float, PI},
    molecule::Molecule,
    vec3d::Vec3d,
};



// ----- PART 1: CGF -----

/// Contracted Gaussian Function (STO-3G),
/// which is linear combination of three primitive gaussian functions.
#[derive(Debug, Clone, Copy)]
pub struct CGF {
    // TODO?: rewrite to `ndarray`
    co: [float; 3],
    alpha: [float; 3],
    coordinate: Vec3d, // TODO: rename to `position`
}
impl CGF {
    pub fn from(zeta: float, n: i32, coordinate: Vec3d) -> Self {
        let constract_cos: [[float; 3]; 2] = [
            [0.444635, 0.535328, 0.154329],
            [0.700115, 0.399513, -0.0999672],
        ];
        let alphas: [[float; 3]; 2] = [
            [0.109818, 0.405771, 2.22766],
            [0.0751386, 0.231031, 0.994203]
        ];
        let n: usize = n as usize;
        let co:    [float; 3] = constract_cos[n-1];
        let alpha: [float; 3] = alphas[n-1].iter()
            .map(|a| a * zeta.powi(2))
            .collect::<Vec<float>>().try_into().unwrap();
        // TODO?: gtos, gto, cgf
        Self {
            co,
            alpha,
            coordinate,
        }
    }

    // fn get_gto() -> ! {
    //     todo!()
    // }
}



// ----- PART 2: Compute integrals between two primitive gaussian functions  -----
// References:
// 1. https://github.com/aced125/Hartree_Fock_jupyter_notebook
// 2. Szabo, Ostlund, Modern Quantum Chemistry, 411-416, https://yyrcd-1256568788.cos.na-siliconvalley.myqcloud.com/yyrcd/2019-11-13-Gaussian_Integral.pdf


// struct GaussFunction {
//     alpha: float,
//     coordinate: float,
// }
// TODO: refactor to `GaussFunction` struct
type GaussFunction = (float, Vec3d);

/// The product of two Gaussians gives another Gaussian. (pp411)
///
/// INPUT:
/// A: (gaussian_A alpha, gaussian_A coordinate)
/// B: (gaussian_B alpha, gaussian_B coordinate)
///
/// OUTPUT:
/// p: New gaussian's alpha
/// diff: squared difference of the two coordinates
/// K: New gaussian's prefactor
/// Rp: New gaussian's coordinate
fn gauss_product(a: GaussFunction, b: GaussFunction) -> (float, float, float, Vec3d) {
    let (a, ra) = a;
    let (b, rb) = b;
    let p = a + b;
    let diff = (ra - rb).len2();
    let n = (4.0*a*b / PI.powi(2)).powf(0.75);
    let k = n * (-a*b/p*diff).exp();
    let rp = (a*ra + b*rb) / p;
    (p, diff, k, rp)
}

/// Compute overlap integral between two primitive gaussian functions.
///
/// INPUT:
/// A: (gaussian_A alpha, gaussian_A coordinate)
/// B: (gaussian_B alpha, gaussian_B coordinate)
fn overlap(a: GaussFunction, b: GaussFunction) -> float {
    let (p, _diff, k, _rp) = gauss_product(a, b);
    let prefactor = (PI / p).powf(1.5);
    prefactor * k
}

/// Compute kinetic integral between two primitive gaussian functions.
///
/// INPUT:
/// A: (gaussian_A alpha, gaussian_A coordinate)
/// B: (gaussian_B alpha, gaussian_B coordinate)
fn kinetic(a: GaussFunction, b: GaussFunction) -> float {
    let (p, diff, k, _rp) = gauss_product(a, b);
    let prefactor = (PI / p).powf(1.5);
    let (a, _ra) = a;
    let (b, _rb) = b;
    let reduced_exponent = a * b / p;
    reduced_exponent * (3.0 - 2.0 * reduced_exponent * diff) * prefactor * k
}

/// Fo function for calculating potential and e-e repulsion integrals.
/// Just a variant of the error function
fn fo(t: float) -> float {
    if t.abs() <= 1e-10 { return 1.0; }
    (0.5 * (PI / t).sqrt()) * erf(t.sqrt())
}

/// Compute Nuclear-electron attraction integral.
///
/// INPUT:
/// A: (gaussian_A alpha, gaussian_A coordinate)
/// B: (gaussian_B alpha, gaussian_B coordinate)
/// coordinate: coordinate of nuclear
/// charge: charge of nuclear
fn potential(a: GaussFunction, b: GaussFunction, coordinate: Vec3d, charge: float) -> float {
    let (p, _diff, k, rp) = gauss_product(a, b);
    let rc = coordinate;
    let zc = charge;
    (-2.0 * PI * zc / p) * k * fo(p * (rp - rc).len2())
}

/// Compute electron-electron repulsion integral.
fn repulsion(
    a: GaussFunction,
    b: GaussFunction,
    c: GaussFunction,
    d: GaussFunction,
) -> float {
    let (p, _diff_ab, k_ab, rp) = gauss_product(a, b);
    let (q, _diff_cd, k_cd, rq) = gauss_product(c, d);
    // TODO: optimize all powf -> powi + sqrt
    let repul_prefactor = 2.0 * PI.powf(2.5) * (p * q * (p + q).sqrt()).powi(-1);
    repul_prefactor*k_ab*k_cd*fo(p*q/(p+q)*(rp-rq).len2())
}



// ----- PART 3: Compute integrals between two contracted gaussian functions -----

/// Compute overlap integral between two contracted gaussian functions.
fn s_int(cgf1: &CGF, cgf2: &CGF) -> float {
    let mut s: float = 0.0;
    // TODO?: rewrite from `.enumerate()` to plain ranges?
    for (i, _) in cgf1.alpha.iter().enumerate() {
        for (j, _) in cgf2.alpha.iter().enumerate() {
            s += cgf1.co[i] * cgf2.co[j] * overlap(
                (cgf1.alpha[i], cgf1.coordinate),
                (cgf2.alpha[j], cgf2.coordinate)
            );
        }
    }
    s
}

/// Compute kinetics integral between two contracted gaussian functions.
fn t_int(cgf1: &CGF, cgf2: &CGF) -> float {
    let mut t: float = 0.0;
    for (i, _) in cgf1.alpha.iter().enumerate() {
        for (j, _) in cgf2.alpha.iter().enumerate() {
            t += cgf1.co[i] * cgf2.co[j] * kinetic(
                (cgf1.alpha[i], cgf1.coordinate),
                (cgf2.alpha[j], cgf2.coordinate)
            );
        }
    }
    t
}


/// Compute electron-nuclear integral between two contracted gaussian functions.
fn v_en_int(cgf1: &CGF, cgf2: &CGF, mol: &Molecule) -> float {
    let mut v: float = 0.0;
    for (i, _) in cgf1.alpha.iter().enumerate() {
        for (j, _) in cgf2.alpha.iter().enumerate() {
            for k in 0..mol.get_number_of_atoms() {
                v += cgf1.co[i] * cgf2.co[j] * potential(
                    (cgf1.alpha[i], cgf1.coordinate),
                    (cgf2.alpha[j], cgf2.coordinate),
                    mol.atoms[k].get_position(),
                    mol.atoms[k].get_charge()
                );
            }
        }
    }
    v
}

/// Compute electron-electron repulsion integral.
fn r_int(cgf1: &CGF, cgf2: &CGF, cgf3: &CGF, cgf4: &CGF) -> float {
    let mut repul: float = 0.0;
    for (i1, _) in cgf1.alpha.iter().enumerate() {
        for (i2, _) in cgf2.alpha.iter().enumerate() {
            for (i3, _) in cgf3.alpha.iter().enumerate() {
                for (i4, _) in cgf4.alpha.iter().enumerate() {
                    let rp = repulsion(
                        (cgf1.alpha[i1], cgf1.coordinate),
                        (cgf2.alpha[i2], cgf2.coordinate),
                        (cgf3.alpha[i3], cgf3.coordinate),
                        (cgf4.alpha[i4], cgf4.coordinate)
                    );
                    repul += cgf1.co[i1] * cgf1.co[i2] * cgf1.co[i3] * cgf1.co[i4] * rp;
                }
            }
        }
    }
    repul
}



// ----- PART 4: Build matrices -----

type Matrix = DMatrix<float>;

fn square_matrix(d: usize) -> Matrix {
    DMatrix::from_element(d, d, 0.0)
}

/// Compute overlap matrix S.
///
/// INPUT:
///     cgfs: basis functions
/// OUTPUT:
///     S: Overlap matrix
fn s_matrix(cgfs: &Vec<CGF>) -> Matrix {
    let mut s = square_matrix(cgfs.len());
    for (i, cgf1) in cgfs.iter().enumerate() {
        for (j, cgf2) in cgfs.iter().enumerate() {
            s[(i, j)] = s_int(cgf1, cgf2);
        }
    }
    s
}


/// Compute the core hamiltonian matrix H.
/// H_core = electron kinetics energy + electron nuclear potential energy
///
/// INPUT:
///     cgfs: basis functions
///     mol: which contain the nuclear charge and nuclear coordinate information
/// OUTPUT:
///     H: core hamiltonian matrix
fn h_matrix(cgfs: &Vec<CGF>, mol: &Molecule) -> Matrix {
    let mut t = square_matrix(cgfs.len());
    let mut v = square_matrix(cgfs.len());
    for (i, cgf1) in cgfs.iter().enumerate() {
        for (j, cgf2) in cgfs.iter().enumerate() {
            t[(i, j)] = t_int(cgf1, cgf2);
            v[(i, j)] = v_en_int(cgf1, cgf2, mol);
        }
    }
    let h = t + v;
    h
}

#[derive(Debug, Clone)]
struct Array4d {
    data: DMatrix<Matrix>
}
impl Array4d {
    pub fn new(d: usize, value: float) -> Self {
        Self {
            data: DMatrix::from_element(
                d, d,
                DMatrix::from_element(
                    d, d,
                    value
                )
            )
        }
    }
    pub fn zeros(d: usize) -> Self { Self::new(d, 0.0) }
    // pub fn shape(&self) -> (usize, usize) { self.data.shape() }
}
impl Index<(usize, usize, usize, usize)> for Array4d {
    type Output = float;
    fn index(&self, (i1, i2, i3, i4): (usize, usize, usize, usize)) -> &Self::Output {
        &self.data[(i1, i2)][(i3, i4)]
    }
}
impl IndexMut<(usize, usize, usize, usize)> for Array4d {
    fn index_mut(&mut self, (i1, i2, i3, i4): (usize, usize, usize, usize)) -> &mut Self::Output {
        &mut self.data[(i1, i2)][(i3, i4)]
    }
}

/// Compute the electron repulsion integral matrix R.
///
/// INPUT:
///     cgfs: basis functions
/// OUTPUT:
///     R: repulsion matrix
fn r_matrix(cgfs: &Vec<CGF>) -> Array4d {
    // start = time.time()
    let mut r = Array4d::zeros(cgfs.len());
    for (i1, cgf1) in cgfs.iter().enumerate() {
        for (i2, cgf2) in cgfs.iter().enumerate() {
            for (i3, cgf3) in cgfs.iter().enumerate() {
                for (i4, cgf4) in cgfs.iter().enumerate() {
                    r[(i1, i2, i3, i4)] = r_int(cgf1, cgf2, cgf3, cgf4);
                }
            }
        }
    }
    // stop = time.time()
    // print('time Repu: {:.1f} s'.format(stop-start))
    r
}

/// Compute density matrix P.
///
/// INPUT:
///     Co: coefficents matrix
///     N: num of electrons
/// OUTPUT:
///     P: repulsion matrix
fn p_matrix(co: &Matrix, n: usize) -> Matrix {
    let mut p = square_matrix(co.shape().0);
    for i1 in 0..co.shape().0 {
        for i2 in 0..co.shape().0 {
            for j in 0..n/2 {
                p[(i1, i2)] += 2.0 * co[(i1, j)]*co[(i2, j)];
            }
        }
    }
    p
}


/// Compute G matrix.
/// G = coulombic repulsion energy + exchange energy
///
/// INPUT:
///     P: density matrix
///     R: electron repulsion matrix
/// OUTPUT:
///     G: repulsion matrix
fn g_matrix(p: &Matrix, r: &Array4d) -> Matrix {
    let num_bfs: usize = p.shape().0;
    let mut g = square_matrix(num_bfs);
    for i1 in 0..num_bfs {
        for i2 in 0..num_bfs {
            let mut g_rs: float = 0.0;
            for i3 in 0..num_bfs {
                for i4 in 0..num_bfs {
                    let int1 = r[(i1, i2, i3, i4)];
                    let int2 = r[(i1, i4, i3, i2)];
                    g_rs += p[(i3, i4)] * (int1 - 0.5 * int2);
                }
            }
            g[(i1, i2)] = g_rs;
        }
    }
    g
}


// /// Compute fock matrix F.
// /// F = H_core + G
// fn f_matrix(h: Matrix, g: Matrix) -> Matrix {
//     h + g
// }



// ----- PART 5: Other Equations -----

/// Compute Nuclear-Nuclear repulsion energy
fn v_nn(mol: &Molecule) -> float {
    let mut nn: float = 0.0;
    for i in 0..mol.get_number_of_atoms() {
        for j in i+1..mol.get_number_of_atoms() {
            // Select atoms from molecule
            let ri = mol.atoms[i].get_position();
            let rj = mol.atoms[j].get_position();
            let zi = mol.atoms[i].get_charge();
            let zj = mol.atoms[j].get_charge();
            nn += zi * zj / (ri - rj).len();
        }
    }
    nn
}

/// Slove secular equation, return the MO energies (eigenvalue) and improved coeffients (eigenvector)
///
/// INPUT:
///     F: fock matrix
///     S: overlap integral
/// OUTPUT:
///     ei: eigenvalue
///     C: eigenvector
fn secular_eqn(f: &Matrix, s: &Matrix) -> (Vec<Complex<float>>, Matrix) {
    todo!()
    // let (ei, c) = eigh(f, s);
    // (ei, c)
}

/// Compute the total energy.
///
/// INPUT:
/// e: MO energies
/// N: num of electrons
/// P: density matrix
/// H: h_core matrix
/// Vnn: nuclear nuclear repulsion energy
fn energy_tot(
    e: &Vec<Complex<float>>,
    n: usize,
    p: &Matrix,
    h: &Matrix,
    vnn: float,
) -> float {
    let mut e_tot: float = 0.0;
    for i in 0..n/2 {
        e_tot += e[i].re;
    }
    // TODO: chech if `(p*h).sum` is correct
    e_tot += 0.5 * (p * h).sum() + vnn;
    e_tot
}



// ----- PART 6: Utils -----

/// Print information while doing SCF interations.
fn print_info(
    s: &Matrix,
    h: &Matrix,
    e: &Vec<Complex<float>>,
    co: &Matrix,
    p: &Matrix,
    hf_e: float,
    // start,
    // stop,
    delta_e: Option<float>,
    verbose: bool,
) {
    let delta_e: float = delta_e.unwrap_or(0.0);
    if verbose {
        println!("Overlap:\n{s}");
        println!("Core hamiltonian:\n{h}");
        println!("Coefficients:\n{co}");
        println!("Density matrix:\n{p}");
        println!("MO energies:");
        let message = e.iter().enumerate()
            .map(|(i, x)| format!("e{} = {:0.3}", i+1, x))
            .collect::<Vec<String>>()
            .join(", ");
        println!("{message}");
    }
    println!("HF energy: {:0.5} (hartree) = {:0.5} (eV)", hf_e, hf_e*27.211);
    if delta_e != 0.0 {
        println!("dE       : {:.2e}", delta_e);
    }
    // println!("time used: {:.1f} s", stop-start);
}


/// Compare calculated result with reference data.
pub fn compare(calculated: float, reference: float, tolerance: Option<float>) {
    let tolerance: float = tolerance.unwrap_or(1.0e-4);
    let delta = (reference - calculated).abs();
    let message = if delta < tolerance { "PASSED" } else { "FAILED" };
    println!("{}{}{}", "-".repeat(32), message, "-".repeat(33));
    println!("cal: {:.7}, ref: {:.7}\n\n", calculated, reference);
}



// ----- PART 7: Final function -----
/// Run restricted hartree fock for 2-electron diatomic molecule.
///
/// INPUT:
///     mol: Molecule
pub fn run_hf(mol: &Molecule) -> float {
    println!("------------------------------ Initialization ------------------------------");
    println!("------------------------- Ignore repulsion integral ------------------------");
    // num of electron
    let n = mol.get_number_of_atoms();

    // initialization
    // start = time.time();
    let h = h_matrix(&mol.cgfs, mol);
    let s = s_matrix(&mol.cgfs);
    let (mut e, mut co) = secular_eqn(&h, &s);
    let mut p = p_matrix(&co, n);
    let vnn = v_nn(mol);
    let mut hf_e = energy_tot(&e, n, &p, &h, vnn);
    // stop = time.time();

    print_info(&s, &h, &e, &co, &p, hf_e, None, false);
    println!("------- Caculating Electron Repulsion Integral (may take some time) --------");
    let r = r_matrix(&mol.cgfs);
    let mut delta_e: float = 1.0;
    let mut iter: u32 = 0;
    let mut previous_e = hf_e;

    // Iterations
    const MAXITER: u32 = 40;
    const E_CONV: float = 1.0e-6;
    while delta_e > E_CONV && iter < MAXITER {
        println!("------------------------------ Iteration {ip1} ------------------------------", ip1=iter+1);
        // start = time.time();

        // important scf steps
        let g = g_matrix(&p, &r);
        let f = &h + g; // TODO: use `f_matrix` function
        (e, co) = secular_eqn(&f, &s);
        p = p_matrix(&co, n);
        hf_e = energy_tot(&e, n, &p, &h, vnn);

        delta_e = (hf_e - previous_e).abs();
        previous_e = hf_e;
        iter += 1;
        // stop = time.time();
        print_info(&s, &h, &e, &co, &p, hf_e, None, false);
    }

    hf_e
}
