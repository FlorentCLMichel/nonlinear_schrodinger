#[allow(unused_imports)]
use crate::{ Complex, FtStruct, MFtStruct, complex::rmse };

#[test]
fn new_struct() {
    let n: usize = 10;
    let _ft_struct = FtStruct::new(n);
}

#[test]
fn ft_1() {
    let n: usize = 10;
    let x = vec![Complex {real: 1., imag: 0.}; n];
    let ft_struct = FtStruct::new(n);
    let y = ft_struct.ft(&x).unwrap();
    let mut y_expected = vec![Complex { real: 0., imag: 0. }; n];
    y_expected[0] = Complex { real: n as f64, imag: 0. };
    assert!(rmse(&y, &y_expected) < 1e-8);
}

#[test]
fn ift_1() {
    let n: usize = 10;
    let x = vec![Complex {real: 1., imag: 0.}; n];
    let ft_struct = FtStruct::new(n);
    let y = ft_struct.ift(&x).unwrap();
    let mut y_expected = vec![Complex { real: 0., imag: 0. }; n];
    y_expected[0] = Complex { real: 1., imag: 0. };
    assert!(rmse(&y, &y_expected) < 1e-8);
}

#[test]
fn ft_2() {
    let n: usize = 10;
    let mut x = vec![Complex {real: 0., imag: 0.}; n];
    x[0] = Complex {real: 1., imag: 0.};
    let ft_struct = FtStruct::new(n);
    let y = ft_struct.ft(&x).unwrap();
    let y_expected = vec![Complex { real: 1., imag: 0. }; n];
    assert!(rmse(&y, &y_expected) < 1e-8);
}

#[test]
fn ift_2() {
    let n: usize = 10;
    let mut x = vec![Complex {real: 0., imag: 0.}; n];
    x[0] = Complex {real: n as f64, imag: 0.};
    let ft_struct = FtStruct::new(n);
    let y = ft_struct.ift(&x).unwrap();
    let y_expected = vec![Complex { real: 1., imag: 0. }; n];
    assert!(rmse(&y, &y_expected) < 1e-8);
}

#[test]
fn ft_3() {
    let n: usize = 4;
    let x = vec![Complex {real: 1., imag: 4.}, 
                 Complex {real: 2., imag: 1.},
                 Complex {real: 3., imag: 4.},
                 Complex {real: 4., imag: 1.}];
    let ft_struct = FtStruct::new(n);
    let y = ft_struct.ft(&x).unwrap();
    let y_expected = vec![Complex { real: 10., imag: 10. },
                          Complex { real: -2., imag: 2. },
                          Complex { real: -2., imag: 6. },
                          Complex { real: -2., imag: -2. }];
    assert!(rmse(&y, &y_expected) < 1e-8);
}

#[test]
fn ift_3() {
    let n: usize = 4;
    let x = vec![Complex { real: 10., imag: 10. },
                 Complex { real: -2., imag: 2. },
                 Complex { real: -2., imag: 6. },
                 Complex { real: -2., imag: -2. }];
    let ft_struct = FtStruct::new(n);
    let y = ft_struct.ift(&x).unwrap();
    let y_expected = vec![Complex {real: 1., imag: 4.}, 
                          Complex {real: 2., imag: 1.},
                          Complex {real: 3., imag: 4.},
                          Complex {real: 4., imag: 1.}];
    assert!(rmse(&y, &y_expected) < 1e-8);
}

#[test]
fn ft_4() {
    let n: usize = 100;
    let mut x = vec![Complex {real: 0., imag: 0.}; n];
    x[0] = Complex {real: 1., imag: 0.};
    let ft_struct = FtStruct::new(n);
    let y = ft_struct.ft(&x).unwrap();
    let y_expected = vec![Complex { real: 1., imag: 0. }; n];
    assert!(rmse(&y, &y_expected) < 1e-8);
}

#[test]
fn ift_4() {
    let n: usize = 100;
    let x = vec![Complex {real: 0., imag: 1.}; n];
    let ft_struct = FtStruct::new(n);
    let y = ft_struct.ift(&x).unwrap();
    let mut y_expected = vec![Complex { real: 0., imag: 0. }; n];
    y_expected[0] = Complex { real: 0., imag: 1. };
    assert!(rmse(&y, &y_expected) < 1e-8);
}

#[test]
fn ft_ift_1() {
    let x = vec![
        Complex { real: 1., imag: 1. },
        Complex { real: 2., imag: 2. },
        Complex { real: 3., imag: 1. },
        Complex { real: 4., imag: 2. },
        Complex { real: 5., imag: 3. },
        Complex { real: 6., imag: 2. },
        Complex { real: 7., imag: 1. },
        Complex { real: 8., imag: 0. },
    ];
    let ft_struct = FtStruct::new(x.len());
    let y = ft_struct.ft(&x).unwrap();
    let z = ft_struct.ift(&y).unwrap();
    assert!(rmse(&x, &z) < 1e-8);
}

#[test]
fn m_ft_1() {
    let shape: Vec<usize> = vec![4,4];
    let x = vec![
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
    ];
    let ft_struct = MFtStruct::new(shape);
    let y = ft_struct.ft(&x).unwrap();
    let y_expected = vec![
        Complex { real: 40., imag: 40. },
        Complex { real: -8., imag: 8. },
        Complex { real: -8., imag: 24. },
        Complex { real: -8., imag: -8. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
    ];
    assert!(rmse(&y, &y_expected) < 1e-8);
}


#[test]
fn m_ft_2() {
    let shape: Vec<usize> = vec![3,4];
    let x = vec![
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
    ];
    let ft_struct = MFtStruct::new(shape);
    let y = ft_struct.ft(&x).unwrap();
    let y_expected = vec![
        Complex { real: 30., imag: 30. },
        Complex { real: -6., imag: 6. },
        Complex { real: -6., imag: 18. },
        Complex { real: -6., imag: -6. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
    ];
    assert!(rmse(&y, &y_expected) < 1e-8);
}


#[test]
fn m_ft_3() {
    let shape: Vec<usize> = vec![5,4];
    let x = vec![
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
    ];
    let ft_struct = MFtStruct::new(shape);
    let y = ft_struct.ft(&x).unwrap();
    let y_expected = vec![
        Complex { real: 50., imag: 50. },
        Complex { real: -10., imag: 10. },
        Complex { real: -10., imag: 30. },
        Complex { real: -10., imag: -10. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
        Complex { real: 0., imag: 0. },
    ];
    assert!(rmse(&y, &y_expected) < 1e-8);
}


#[test]
fn m_ft_ift_1() {
    let shape: Vec<usize> = vec![4,4];
    let x = vec![
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
    ];
    let ft_struct = MFtStruct::new(shape);
    let y = ft_struct.ft(&x).unwrap();
    let y = ft_struct.ift(&y).unwrap();
    assert!(rmse(&y, &x) < 1e-8);
}


#[test]
fn m_ft_ift_2() {
    let shape: Vec<usize> = vec![2,4,2];
    let x = vec![
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
        Complex {real: 1., imag: 4.}, 
        Complex {real: 2., imag: 1.},
        Complex {real: 3., imag: 4.},
        Complex {real: 4., imag: 1.},
    ];
    let ft_struct = MFtStruct::new(shape);
    let y = ft_struct.ft(&x).unwrap();
    let y = ft_struct.ift(&y).unwrap();
    assert!(rmse(&y, &x) < 1e-8);
}
