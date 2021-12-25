use crate::{ R, C, PI };

#[cfg(feature="multithread_ft")]
use rayon::prelude::*;


/// find the optimal plan for a given length
///
/// # Argument
///
/// * `n`: the Fourier transform size
///
/// # Return 
///
/// A vector of `usize`s with the size of the Fourier transform to perform at each step.
///
/// # Example
///
/// ```
/// use nonlinear_schrodinger::*;
///
/// # fn main() {
/// # 
/// let n: usize = 100;
/// let plan = find_plan(n);
///
/// assert_eq!(plan, vec![2, 2, 5, 5]);
/// # 
/// # }
/// ```
pub fn find_plan(n: usize) -> Vec<usize> {
    let mut plan = Vec::<usize>::new();
    let mut n = n;
    let mut m: usize = 2;
    let mut step: usize = 1;
    while n > 1 {
        while n%m == 0 {
            plan.push(m);
            n /= m;
        }
        m += step;
        if step==1 {
            step = 2;
        }
    }
    plan
}


/// get the butterfly coefficients 
///
/// # Arguments
///
/// * `n`: the Fourier transform size
/// * `plan`: the Fourier transform plan (*e.g.* as returned by [`find_plan`])
///
/// # Return 
///
/// A triplet of vectors of the form `(v_butterfly, v_twiddle_global, v_twiddle_local)` where
///
/// * `v_butterfly` contains the butterfly coefficients,
/// * `v_twiddle_global` contains the twiddle factors applied during the butterfly,
/// * `v_twiddle_local` contains the twiddle factors for the smaller Fourier transforms.
///
/// # Example
///
/// ```
/// use nonlinear_schrodinger::*;
///
/// # fn main() {
/// # 
/// let n: usize = 100;
/// let plan = find_plan(n);
/// let (butterfly, twiddle_global, twidle_local) = get_butterflies_and_twiddles(n, &plan);
/// # 
/// # }
/// ```
pub fn get_butterflies_and_twiddles(n: usize, plan: &[usize]) 
    -> (Vec<usize>, Vec<C>, Vec<C>)
{
    let mut butterflies: Vec<Vec<usize>> = vec![(0..n).collect()];
    let mut twiddles: Vec<Vec<C>> = vec![vec![C::new(1.,0.); n]];
    let mut twiddles_small_fts = Vec::<Vec<C>>::new();
    get_butterflies_and_twiddles_(n, plan, 1, &mut butterflies, &mut twiddles, &mut twiddles_small_fts);
    (butterflies.into_iter().flatten().collect(), 
     twiddles.into_iter().flatten().collect(), 
     twiddles_small_fts.into_iter().flatten().collect())
}


/// A naive Fourier transform function with pre-computed twiddle factors
/// 
/// # Arguments 
///
/// * `a`: input vector
/// * `b`: output vector
/// * `n`: size of the Fourier transforms
/// * `tf`: twiddle factors
///
/// The result is stored in the array `b`.
/// The vectors `a` and `b` must have the same length, which must be a multiple of `n`.
///
/// # Example
///
/// ```
/// use nonlinear_schrodinger::*;
///
/// # fn main() {
/// # 
/// // size of the Fourier transform
/// let n: usize = 4;
///
/// // twiddle factors
/// let tf = vec![C::new(1.,0.), C::new(0.,1.), C::new(-1.,0.), C::new(0.,-1.)];
///
/// // input vector
/// let x = vec![C::new(1., 2.); 4];
///
/// // vector to store the output
/// let mut y = vec![C::new(0., 0.); 4];
///
/// // Fourier transform
/// ft_inplace_tf(&x, &mut y, n, &tf);
///
/// // check the result
/// let expected_y= vec![C::new(4.,8.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.)];
/// assert_eq!(y, expected_y);
/// # 
/// # }
/// ```
pub fn ft_inplace_tf (a: &[C], b: &mut [C], n: usize, tf: &[C]) {
    
    #[cfg(feature="multithread_ft")]
    {
        a.par_iter().chunks(n).zip(b.par_iter_mut().chunks(n)).for_each(move |(a_c, mut b_c)| {
            let mut coeff: C;
            let mut index_tf: usize;
            for (i, x_b) in b_c.iter_mut().enumerate().take(n) {
                coeff = C::new(0.,0.);
                index_tf = 0;
                for &x_a in a_c.iter() {
                    coeff += tf[index_tf] * *x_a;
                    index_tf += i;
                    if index_tf >= n { index_tf -= n; }
                }
                **x_b = coeff;
            }
        });
    } 
    
    #[cfg(not(feature="multithread_ft"))]
    {
        a.chunks(n).zip(b.chunks_mut(n)).for_each(|(a_c, b_c)| {
            let mut coeff: C;
            let mut index_tf: usize;
            for (i, x_b) in b_c.iter_mut().enumerate().take(n) {
                coeff = C::new(0.,0.);
                index_tf = 0;
                for &x_a in a_c {
                    coeff += tf[index_tf] * x_a;
                    index_tf += i;
                    if index_tf >= n { index_tf -= n; }
                }
                *x_b = coeff;
            }
        });
    }
}


pub fn ft_inplace_pow2 (a: &[C], b: &mut [C]) {
    
    #[cfg(feature="multithread_ft")]
    {
        a.par_iter().chunks(2).zip(b.par_iter_mut().chunks(2)).for_each(|(a_c, mut b_c)| {
            *b_c[0] = *a_c[0] + *a_c[1];
            *b_c[1] = *a_c[0] - *a_c[1];
        });
    }
    
    #[cfg(not(feature="multithread_ft"))]
    {
        a.chunks(2).zip(b.chunks_mut(2)).for_each(|(a_c, b_c)| {
            b_c[0] = a_c[0] + a_c[1];
            b_c[1] = a_c[0] - a_c[1];
        });
    }
}


fn get_butterflies_and_twiddles_(n: usize, plan: &[usize], g: usize, 
                                 butterflies: &mut Vec<Vec<usize>>,
                                 twiddles: &mut Vec<Vec<C>>,
                                 twiddles_small_ft: &mut Vec<Vec<C>>) 
{
   
    // if the plan is empty, do nothing
    if plan.is_empty() {
        return;
    }
    
    // if the plan has a length 1, perform the last operations and stop
    if plan.len() == 1 {
        butterflies.push((0..n).collect());
        twiddles.push(vec![C::new(1.,0.); n]);
        twiddles_small_ft.push(compute_twiddles(plan[0]));
        return;
    }

    // number of vectors already in the array
    let mut last_index = twiddles.len() - 1;

    // divide a into sub-arrays
    let previous_butterfly = butterflies[last_index].clone();
    let d = plan[0];
    let m = n / (d*g);
    for k in 0..g {
        for i in 0..d {
            for j in 0..m {
                butterflies[last_index][j+i*m+k*d*m] = previous_butterfly[i+j*d+k*d*m];
            }
        }
    }

    get_butterflies_and_twiddles_(n, &plan[1..], g*d, butterflies, twiddles, twiddles_small_ft);
    twiddles_small_ft.push(compute_twiddles(plan[0]));
    
    // update the number of vectors already in the array
    last_index = twiddles.len() - 1;

    // vector of complex exponentials to use
    let mut exps = Vec::<C>::with_capacity(d*m);
    let d_m_r = (d*m) as R;
    for i in 0..(d*m) {
        exps.push(C::new(0.,-2.*PI*(i as R) / d_m_r).exp());
    }

    // transposition and twiddle factors
    let previous_butterfly = butterflies[last_index].clone();
    for k in 0..g {
        let k_d_m = k*d*m;
        for i in 0..m {
            for j in 0..d {
                butterflies[last_index][j+i*d+k_d_m] = previous_butterfly[i+j*m+k_d_m];
                twiddles[last_index][j+i*d+k_d_m] *= exps[(i*j) % (m*d)];
            }
        }
    }

    // second Fourier transform
    butterflies.push((0..n).collect());
    twiddles.push(vec![C::new(1.,0.); n]);
    last_index += 1;

    // transpose
    let previous_butterfly = butterflies[last_index].clone();
    for k in 0..g {
        let k_d_m = k*d*m;
        for i in 0..d {
            for j in 0..m {
                butterflies[last_index][j+i*m+k_d_m] = previous_butterfly[i+j*d+k_d_m];
            }
        }
    }
}


fn compute_twiddles(n: usize) -> Vec<C> {
    let mut coeffs = Vec::<C>::with_capacity(n);
    let x: R = -2.*PI / (n as R);
    let c = C::new(0.,x);
    for i in 0..n {
        coeffs.push((c*(i as R)).exp());
    }
    coeffs
}
