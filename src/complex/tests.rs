use crate::Complex;

#[test]
fn add_1() {
    let z1 = Complex {real: 0., imag: 1.};
    let z2 = Complex {real: 1., imag: 0.};
    let z3 = Complex {real: 1., imag: 1.};
    assert_eq!(z1 + z2, z3);
}

#[test]
fn add_2() {
    let mut z1 = Complex {real: 0., imag: 1.};
    z1 += Complex {real: 1., imag: -1.};
    let z2 = Complex {real: 1., imag: 0.};
    assert_eq!(z1, z2);
}

#[test]
fn mul_1() {
    let z1 = Complex {real: 0., imag: 1.};
    let z2 = Complex {real: 1., imag: -1.};
    let z3 = Complex {real: 1., imag: 1.};
    assert_eq!(z1 * z2, z3);
}

#[test]
fn mul_2() {
    let mut z1 = Complex {real: 0., imag: 1.};
    z1 *= Complex {real: 1., imag: -1.};
    let z2 = Complex {real: 1., imag: 1.};
    assert_eq!(z1, z2);
}

#[test]
fn mul_3() {
    let mut z1 = Complex {real: -2., imag: 1.};
    z1 *= Complex {real: 1., imag: -2.};
    let z2 = Complex {real: 0., imag: 5.};
    assert_eq!(z1, z2);
}

#[test]
fn exp_1() {
    let z1 = Complex {real: 0., imag: 0.};
    let z2 = z1.exp();
    let z3 = Complex {real: 1., imag: 0.};
    assert_eq!(z2, z3);
}

#[test]
fn exp_2() {
    let z1 = Complex {real: 1., imag: 1.};
    let z2 = z1.exp();
    let z3 = Complex {real: 1_f64.exp()*1_f64.cos(), imag: 1_f64.exp()*1_f64.sin()};
    assert_eq!(z2, z3);
}

#[test]
fn abs_1() {
    let z1 = Complex {real: -4., imag: 3.};
    let r1 = 5.;
    assert_eq!(z1.abs(), r1);
}

#[test]
fn abs_2() {
    let z1 = Complex {real: 1., imag: 1.};
    let r1 = 2_f64.sqrt();
    assert_eq!(z1.abs(), r1);
}
