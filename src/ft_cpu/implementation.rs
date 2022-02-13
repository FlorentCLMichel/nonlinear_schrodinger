use super::{ R, C, FtStruct, MFtStruct, FFTError };
use crate::PI;
use super::plan::*;
use std::cell::RefCell;

impl FtStruct {

    /// create a new ft structure
    ///
    /// # Argument
    ///
    /// * `n`: size of the Fourier transform
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::*;
    ///
    /// # fn main() {
    /// let ft_size: usize = 10;
    /// let ft_struct = FtStruct::new(ft_size);
    /// # }
    /// ```
    pub fn new(n: usize) -> FtStruct {

        // see if the size is a power of 2
        let is_power_2 = check_pow_2(n);

        // compute c_bar and f_c_tilde if necessary
        let mut c_bar = Vec::<C>::new();
        let mut f_c_tilde = Vec::<C>::new();
        let mut buffer1 = Vec::<C>::new();
        let mut buffer2 = Vec::<C>::new();
        let butterflies: Vec<usize>;
        let twiddles: Vec<C>;
        if is_power_2 {
            let plan = find_plan(n);
            let b_tw = get_butterflies_and_twiddles(n, &plan);
            butterflies = b_tw.0;
            twiddles = b_tw.1;
        } else {
            let mut m = 1;
            while 2*n > m+1 {
                m *= 2;
            }
            let plan = find_plan(m);
            let b_tw = get_butterflies_and_twiddles(m, &plan);
            butterflies = b_tw.0;
            twiddles = b_tw.1;
            c_bar = compute_c_bar(n);
            f_c_tilde = compute_f_c_tilde(n, m);
            buffer1 = vec![C::new(0.,0.); m];
            buffer2 = vec![C::new(0.,0.); m];
        }

        // return the ft structure
        FtStruct { n, butterflies, twiddles, c_bar, f_c_tilde, 
            buffer1: RefCell::new(buffer1), buffer2: RefCell::new(buffer2), 
            is_power_2 }
    }


    /// direct Fourier transform
    ///
    /// # Argument
    ///
    /// * `a`: input vector
    ///
    /// # Return
    ///
    /// `Ok(Vec<C>)` containing the result if the computation is successful or `Err(FFTError)` 
    /// otherwise
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::*;
    /// use nonlinear_schrodinger::complex::rmse;
    ///
    /// # fn main() -> Result<(), FFTError> {
    /// // Fourier transform size
    /// let ft_size: usize = 10;
    ///
    /// // build the Fourier structure
    /// let ft_struct = FtStruct::new(ft_size);
    ///
    /// // input
    /// let x: Vec<C> = vec![C::new(1.,0.); ft_size];
    ///
    /// // Fourier transform
    /// let y = ft_struct.ft(&x)?;
    ///
    /// // check the result
    /// let mut expected_y = vec![C::new(0.,0.); ft_size];
    /// expected_y[0] = C::new(ft_size as R, 0.);
    /// assert!(rmse(&y, &expected_y) < 1e-8);
    /// # Ok(())
    /// # }
    /// ```
    pub fn ft(&self, a: &[C]) 
        -> Result<Vec<C>, FFTError>
    {

        let mut b = a.to_vec();
        let mut c = vec![C::new(0.,0.); a.len()];
        self.ft_inplace(&mut b, &mut c)?;
        Ok(c)
    }

    
    /// inverse Fourier transform
    ///
    /// # Argument
    ///
    /// * `a`: input vector
    ///
    /// # Return
    ///
    /// `Ok(Vec<C>)` containing the result if the computation is successful or `Err(FFTError)` 
    /// otherwise
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::*;
    /// use nonlinear_schrodinger::complex::rmse;
    ///
    /// # fn main() -> Result<(), FFTError> {
    /// // Fourier transform size
    /// let ft_size: usize = 10;
    ///
    /// // build the Fourier structure
    /// let ft_struct = FtStruct::new(ft_size);
    ///
    /// // input
    /// let x: Vec<C> = vec![C::new(1.,0.); ft_size];
    ///
    /// // inverse Fourier transform
    /// let y = ft_struct.ift(&x)?;
    ///
    /// // check the result
    /// let mut expected_y = vec![C::new(0.,0.); ft_size];
    /// expected_y[0] = C::new(1., 0.);
    /// assert!(rmse(&y, &expected_y) < 1e-8);
    /// # Ok(())
    /// # }
    /// ```
    pub fn ift(&self, a: &[C]) 
        -> Result<Vec<C>, FFTError>
    {

        let mut b = a.to_vec();
        let mut c = vec![C::new(0.,0.); a.len()];
        self.ift_inplace(&mut b, &mut c)?;
        Ok(c)
    }


    /// direct inplace Fourier transform
    ///
    /// # Argument
    ///
    /// * `a`: input vector
    /// * `b`: output vector
    ///
    /// # Return
    ///
    /// `Ok(())` if the computation is successful or `Err(FFTError)` otherwise
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::*;
    /// use nonlinear_schrodinger::complex::rmse;
    ///
    /// # fn main() -> Result<(), FFTError> {
    /// // Fourier transform size
    /// let ft_size: usize = 10;
    ///
    /// // build the Fourier structure
    /// let ft_struct = FtStruct::new(ft_size);
    ///
    /// // input
    /// let mut x: Vec<C> = vec![C::new(1.,0.); ft_size];
    ///
    /// // vector to store the output
    /// let mut y: Vec<C> = vec![C::new(1.,0.); ft_size];
    ///
    /// // Fourier transform
    /// ft_struct.ft_inplace(&mut x, &mut y)?;
    ///
    /// // check the result
    /// let mut expected_y = vec![C::new(0.,0.); ft_size];
    /// expected_y[0] = C::new(ft_size as R, 0.);
    /// assert!(rmse(&y, &expected_y) < 1e-8);
    /// # Ok(())
    /// # }
    /// ```
    pub fn ft_inplace (&self, a: &mut [C], b: &mut [C])
        -> Result<(),FFTError> 
    {

        // check that a and b have the right length
        if a.len() != self.n {
            return Err(FFTError::new(format!("Wrong first array size: expected {}, got {}", 
                                             self.n, a.len())));
        }
        if b.len() != self.n {
            return Err(FFTError::new(format!("Wrong second array size: expected {}, got {}", 
                                             self.n, b.len())));
        }
        
        // small Fourier transform
        if self.is_power_2 {
            self.compute_fft_pow2(a, b)
        } else {
            self.compute_fft_nonpow2(a, b)
        }
    }
    

    fn compute_fft_nonpow2 (&self, a: &mut [C], b: &mut [C])
        -> Result<(),FFTError> 
    {

        let mut buffer1 = self.buffer1.borrow_mut();
        let mut buffer2 = self.buffer2.borrow_mut();

        // multiply by c_bar
        for (&x,(&y,z)) in a.iter().zip(self.c_bar.iter().zip(buffer1.iter_mut())) {
            *z = x*y;
        } 

        // fill the rest of the buffer with 0s
        for z in buffer1[a.len()..].iter_mut() {
            *z = C::new(0.,0.);
        }

        // Fourier transform
        self.compute_fft_pow2(&mut buffer1, &mut buffer2)?;

        // multiplication by f_c_tilde
        for (x, &y) in buffer2.iter_mut().zip(self.f_c_tilde.iter()) {
            *x *= y;
        }

        // inverse Fourier transform
        self.compute_fft_pow2(&mut buffer2, &mut buffer1)?;
        ft_to_ift(&mut buffer1);

        // save the result in b
        for (x,(&y,&z)) in b.iter_mut().zip(buffer1.iter().zip(self.c_bar.iter())) {
            *x = y * z;
        }

        Ok(()) 
    }

    fn compute_fft_pow2 (&self, a: &mut [C], b: &mut [C])
        -> Result<(),FFTError> 
    {
    
        let n = a.len();
       
        // first butterfly and twiddle
        b.iter_mut().zip(self.butterflies.iter()).for_each(|(e, &butterfly)| { *e = a[butterfly]; });
    
        for (butterflies, twiddles) in self.butterflies[n..].chunks_exact(n)
                                                       .zip(self.twiddles[n..].chunks_exact(n)) {
            
            // small Fourier transform
            ft_inplace_pow2(b, a);
    
            // butterfly and twiddle
            b.chunks_exact_mut(2).zip(butterflies.chunks_exact(2).zip(twiddles.chunks_exact(2)))
                                 .for_each(
                |(e, (butterfly, twiddle))| { 
                    e[0] = a[butterfly[0]]; 
                    e[1] = a[butterfly[1]] * twiddle[1]; 
                });
        }
            
        Ok(()) 
    }


    /// inverse inplace Fourier transform
    ///
    /// # Argument
    ///
    /// * `a`: input vector
    /// * `b`: output vector
    ///
    /// # Return
    ///
    /// `Ok(())` if the computation is successful or `Err(FFTError)` otherwise
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::*;
    /// use nonlinear_schrodinger::complex::rmse;
    ///
    /// # fn main() -> Result<(), FFTError> {
    /// // Fourier transform size
    /// let ft_size: usize = 10;
    ///
    /// // build the Fourier structure
    /// let ft_struct = FtStruct::new(ft_size);
    ///
    /// // input
    /// let mut x: Vec<C> = vec![C::new(1.,0.); ft_size];
    ///
    /// // vector to store the output
    /// let mut y: Vec<C> = vec![C::new(1.,0.); ft_size];
    ///
    /// // inverse Fourier transform
    /// ft_struct.ift_inplace(&mut x, &mut y)?;
    ///
    /// // check the result
    /// let mut expected_y = vec![C::new(0.,0.); ft_size];
    /// expected_y[0] = C::new(1., 0.);
    /// assert!(rmse(&y, &expected_y) < 1e-8);
    /// # Ok(())
    /// # }
    /// ```
    pub fn ift_inplace (&self, a: &mut [C], b: &mut [C])
        -> Result<(),FFTError> 
    {
        
        let n_r = self.n as R;
        match self.ft_inplace(a, b) {
            Ok(()) => (),
            err => {return err;}
        };
    
        b[0] /= n_r;
        for i in 1..=((self.n-1)/2) {
            let buffer = b[i] / n_r;
            b[i] = b[self.n-i] / n_r;
            b[self.n-i] = buffer;
        }
        if self.n%2 == 0 {
            b[self.n/2] /= n_r;
        }
        
        Ok(())
    }
}


impl MFtStruct {

    /// create a new multi-dimensional ft structure
    ///
    /// # Argument
    ///
    /// * `shape`: shape of the Fourier transform
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::*;
    ///
    /// # fn main() {
    /// let ft_shape: Vec<usize> = vec![2, 3, 4];
    /// let ft_struct = MFtStruct::new(ft_shape);
    /// # }
    /// ```
    pub fn new(shape: Vec<usize>) -> MFtStruct {

        // number of dimensions
        let dimension = shape.len();

        // compute the total length
        let mut total_length: usize = 1;
        for l in shape.iter() {
            total_length *= *l;
        }

        // build te 1D FT structures
        let mut ft_structs = Vec::<FtStruct>::new();
        for &x in &shape {
            ft_structs.push(FtStruct::new(x));
        }

        // return the ft structure
        MFtStruct { dimension, total_length, shape, ft_structs }
    }
    

    /// Fourier transform
    ///
    /// # Argument
    ///
    /// * `x`: array to be FTed
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::*;
    /// use nonlinear_schrodinger::complex::rmse;
    ///
    /// # fn main() -> Result<(), FFTError> {
    /// // create the FT structure
    /// let ft_shape: Vec<usize> = vec![1, 3, 4];
    /// let ft_struct = MFtStruct::new(ft_shape);
    ///
    /// // input
    /// let x = vec![
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    /// ];
    /// 
    /// // perform the Fourier transform
    /// let y = ft_struct.ft(&x)?;
    ///
    /// // check the result
    /// let expected_y = vec![
    ///     C::new(3., 0.), C::new(3.,0.), C::new(3.,0.), C::new(3.,0.),
    ///     C::new(0., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(0., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    /// ];
    /// assert!(rmse(&y, &expected_y) < 1e-8);
    /// 
    /// # Ok(())
    /// # }
    /// ```
    pub fn ft(&self, x: &[C]) -> Result<Vec<C>, FFTError> {
        let mut z = x.to_owned();
        let mut y = vec![C::new(0.,0.); x.len()];
        self.ft_inplace(&mut z, &mut y)?;
        Ok(y)
    }
    

    /// inverse Fourier transform
    ///
    /// # Argument
    ///
    /// * `x`: array to be FTed
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::*;
    /// use nonlinear_schrodinger::complex::rmse;
    ///
    /// # fn main() -> Result<(), FFTError> {
    /// // create the FT structure
    /// let ft_shape: Vec<usize> = vec![1, 3, 4];
    /// let ft_struct = MFtStruct::new(ft_shape);
    ///
    /// // input
    /// let x = vec![
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    /// ];
    /// 
    /// // perform the Fourier transform
    /// let y = ft_struct.ift(&x)?;
    ///
    /// // check the result
    /// let expected_y = vec![
    ///     C::new(0.25, 0.), C::new(0.25,0.), C::new(0.25,0.), C::new(0.25,0.),
    ///     C::new(0., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(0., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    /// ];
    /// assert!(rmse(&y, &expected_y) < 1e-8);
    /// 
    /// # Ok(())
    /// # }
    /// ```
    pub fn ift(&self, x: &[C]) -> Result<Vec<C>, FFTError> {
        let mut z = x.to_owned();
        let mut y = vec![C::new(0.,0.); x.len()];
        self.ift_inplace(&mut z, &mut y)?;
        Ok(y)
    }
    
    
    /// inplace Fourier transform
    ///
    /// # Argument
    ///
    /// * `x`: array to be FTed
    /// * `y`: array to store the result
    ///
    /// # Examples
    ///
    /// ```
    /// use nonlinear_schrodinger::*;
    /// use nonlinear_schrodinger::complex::rmse;
    ///
    /// # fn main() -> Result<(), FFTError> {
    /// // create the FT structure
    /// let ft_shape: Vec<usize> = vec![3, 4];
    /// let ft_struct = MFtStruct::new(ft_shape);
    ///
    /// // input
    /// let mut x = vec![
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    /// ];
    /// 
    /// // perform the Fourier transform
    /// let mut y = vec![C::new(0.,0.); x.len()];
    /// ft_struct.ft_inplace(&mut x, &mut y)?;
    ///
    /// // check the result
    /// let expected_y = vec![
    ///     C::new(3., 0.), C::new(3.,0.), C::new(3.,0.), C::new(3.,0.),
    ///     C::new(0., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(0., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    /// ];
    /// assert!(rmse(&y, &expected_y) < 1e-8);
    /// 
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// ```
    /// use nonlinear_schrodinger::*;
    /// use nonlinear_schrodinger::complex::rmse;
    ///
    /// # fn main() -> Result<(), FFTError> {
    /// // create the FT structure
    /// let ft_shape: Vec<usize> = vec![1, 3, 4];
    /// let ft_struct = MFtStruct::new(ft_shape);
    ///
    /// // input
    /// let mut x = vec![
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    /// ];
    /// 
    /// // perform the Fourier transform
    /// let mut y = vec![C::new(0.,0.); x.len()];
    /// ft_struct.ft_inplace(&mut x, &mut y)?;
    ///
    /// // check the result
    /// let expected_y = vec![
    ///     C::new(3., 0.), C::new(3.,0.), C::new(3.,0.), C::new(3.,0.),
    ///     C::new(0., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(0., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    /// ];
    /// assert!(rmse(&y, &expected_y) < 1e-8);
    /// 
    /// # Ok(())
    /// # }
    /// ```
    #[allow(clippy::many_single_char_names)]
    pub fn ft_inplace(&self, x: &mut [C], y: &mut [C]) -> Result<(), FFTError> {
        
        if self.dimension == 0 {
            return Ok(());
        }

        // last Fourier transform
        let n = self.dimension-1;
        for k in 0..(self.total_length / self.shape[n]) {
            self.ft_structs[n].ft_inplace(&mut x[k*self.shape[n]..(k+1)*self.shape[n]], 
                                          &mut y[k*self.shape[n]..(k+1)*self.shape[n]])?;
        }

        // other Fourier transforms
        let mut a: &mut [C];
        let mut b: &mut [C];
        for i in 0..n {
            if i%2 == 1 {
                a = x;
                b = y;
            } else {
                a = y;
                b = x;
            }
            self.transpose(a, b, i, n, self.shape[i], self.shape[n])?;
            for k in 0..(self.total_length / self.shape[i]) {
                self.ft_structs[i].ft_inplace(&mut b[k*self.shape[i]..(k+1)*self.shape[i]],
                                              &mut a[k*self.shape[i]..(k+1)*self.shape[i]])?;
            }
            self.transpose(a, b, i, n, self.shape[n], self.shape[i])?;
        }

        if self.dimension % 2 == 0 {
            y.clone_from_slice(&x[..]);
        }

        Ok(())
    }
    

    /// inplace inverse Fourier transform
    ///
    /// # Argument
    ///
    /// * `x`: array to be FTed
    /// * `y`: array to store the result
    ///
    /// # Examples
    ///
    /// ```
    /// use nonlinear_schrodinger::*;
    /// use nonlinear_schrodinger::complex::rmse;
    ///
    /// # fn main() -> Result<(), FFTError> {
    /// // create the FT structure
    /// let ft_shape: Vec<usize> = vec![3, 4];
    /// let ft_struct = MFtStruct::new(ft_shape);
    ///
    /// // input
    /// let mut x = vec![
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    /// ];
    /// 
    /// // perform the Fourier transform
    /// let mut y = vec![C::new(0.,0.); x.len()];
    /// ft_struct.ift_inplace(&mut x, &mut y)?;
    ///
    /// // check the result
    /// let expected_y = vec![
    ///     C::new(0.25, 0.), C::new(0.25,0.), C::new(0.25,0.), C::new(0.25,0.),
    ///     C::new(0., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(0., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    /// ];
    /// assert!(rmse(&y, &expected_y) < 1e-8);
    /// 
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// ```
    /// use nonlinear_schrodinger::*;
    /// use nonlinear_schrodinger::complex::rmse;
    ///
    /// # fn main() -> Result<(), FFTError> {
    /// // create the FT structure
    /// let ft_shape: Vec<usize> = vec![3, 1, 4];
    /// let ft_struct = MFtStruct::new(ft_shape);
    ///
    /// // input
    /// let mut x = vec![
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(1., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    /// ];
    /// 
    /// // perform the Fourier transform
    /// let mut y = vec![C::new(0.,0.); x.len()];
    /// ft_struct.ift_inplace(&mut x, &mut y)?;
    ///
    /// // check the result
    /// let expected_y = vec![
    ///     C::new(0.25, 0.), C::new(0.25,0.), C::new(0.25,0.), C::new(0.25,0.),
    ///     C::new(0., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    ///     C::new(0., 0.), C::new(0.,0.), C::new(0.,0.), C::new(0.,0.),
    /// ];
    /// assert!(rmse(&y, &expected_y) < 1e-8);
    /// 
    /// # Ok(())
    /// # }
    /// ```
    #[allow(clippy::many_single_char_names)]
    pub fn ift_inplace(&self, x: &mut [C], y: &mut [C]) -> Result<(), FFTError> {
        
        if self.dimension == 0 {
            return Ok(());
        }

        // last Fourier transform
        let n = self.dimension-1;
        for k in 0..(self.total_length / self.shape[n]) {
            self.ft_structs[n].ift_inplace(&mut x[k*self.shape[n]..(k+1)*self.shape[n]], 
                                           &mut y[k*self.shape[n]..(k+1)*self.shape[n]])?;
        }

        // other Fourier transforms
        let mut a: &mut [C];
        let mut b: &mut [C];
        for i in 0..n {
            if i%2 == 1 {
                a = x;
                b = y;
            } else {
                a = y;
                b = x;
            }
            self.transpose(a, b, i, n, self.shape[i], self.shape[n])?;
            for k in 0..(self.total_length / self.shape[i]) {
                self.ft_structs[i].ift_inplace(&mut b[k*self.shape[i]..(k+1)*self.shape[i]],
                                               &mut a[k*self.shape[i]..(k+1)*self.shape[i]])?;
            }
            self.transpose(a, b, i, n, self.shape[n], self.shape[i])?;
        }

        if self.dimension % 2 == 0 {
            y.clone_from_slice(&x[..]);
        }

        Ok(())
    }


    fn transpose(&self, x: &[C], y: &mut[C], mut i: usize, mut j: usize, l_i:  usize, l_j: usize) 
        -> Result<(), FFTError> 
    {

        // check that x has the right length
        if x.len() != self.total_length {
            return Err(FFTError::new(format!("Invalid input length: expected {}, got {}", 
                                             self.total_length, x.len())));
        }

        // check that i and j are smaller than `self.dimension`
        if (i >= self.dimension) || (j >= self.dimension) {
            return Err(FFTError::new(format!("The row a nd column indices ({},{}) must me smaller than the number of dimensions ({})", 
                                             i, j, self.dimension)));
        }

        // exchange i and j if i > j
        if i > j { std::mem::swap(&mut i, &mut j); }

        // compute the lengths
        let mut product_length_before: usize = 1;
        let mut product_length_inter: usize = 1;
        for k in 0..i {
            product_length_before *= self.shape[k];
        }
        for k in (i+1)..j {
            product_length_inter *= self.shape[k];
        }
        let product_length_before = product_length_before;
        let product_length_inter = product_length_inter;
        let product_length_after = self.total_length / (product_length_before * product_length_inter * l_i * l_j);
        
        // perform the transposition
        for k1 in 0..product_length_before {
            for k2 in 0..l_j {
                for k3 in 0..product_length_inter {
                    for k4 in 0..l_i {
                        for k5 in 0..product_length_after {
                            y[
                                k1*l_j*product_length_inter*l_i*product_length_after
                                + k2*product_length_inter*l_i*product_length_after
                                + k3*l_i*product_length_after
                                + k4*product_length_after
                                + k5
                            ] =
                            x[
                                k1*l_i*product_length_inter*l_j*product_length_after
                                + k4*product_length_inter*l_j*product_length_after
                                + k3*l_j*product_length_after
                                + k2*product_length_after
                                + k5
                            ]
                        }
                    }
                }
            }
        }
        
        Ok(())
    }
}


fn check_pow_2(mut n: usize) -> bool {
    
    if n==0 { return false; }
    
    while n > 1 {
        if n%2 == 1 { return false; }
        n /= 2;
    }
    
    true
}


fn compute_c_bar(size: usize) -> Vec<C> {
    (0..size).map(|i| C::new(0., -((i*i) as R)*PI/(size as R)).exp())
             .collect()
}


fn compute_f_c_tilde(size: usize, size_large_ft: usize) -> Vec<C> {
    
    // compute the vector c_tilde
    let mut c_tilde = vec![C { real: 0., imag: 0. }; size_large_ft];

    #[allow(clippy::needless_range_loop)]
    for i in 0..size {
        c_tilde[i] = C::new(0., ((i*i) as R)*PI/(size as R)).exp();
    }
    
    for i in 1..size {
        c_tilde[size_large_ft-i] = C::new(0., ((i*i) as R)*PI/(size as R)).exp();
    }
    
    // compute the Fourier transform of c_tilde
    let mut f_c_tilde = c_tilde.clone();
    fft_pow2(&mut c_tilde.clone(), &mut f_c_tilde).unwrap();

    f_c_tilde
}


fn fft_pow2(a: &mut [C], b: &mut [C]) 
    -> Result<(),FFTError> 
{
    if a.len() != b.len() {
        return Err(FFTError::new(format!(
            "The arrays a (length {}) and b (length {}) must have the same length",
            a.len(), b.len()
        )));
    }
    if (a.is_empty()) || ((a.len() & (a.len() - 1)) != 0) {
        return Err(FFTError::new(format!(
            "The arrays a (length {}) needs to have a power-of-two length",
            a.len()
        )));
    }
    let n = a.len();
    let plan = find_plan(n);
    let butterflies_and_twiddles = get_butterflies_and_twiddles(n, &plan);
    fft_pow2_from_tb(a, b, &butterflies_and_twiddles.0, &butterflies_and_twiddles.1)
}


fn fft_pow2_from_tb (a: &mut [C], b: &mut [C], 
                     butterflies: &[usize], 
                     twiddles: &[C]) 
    -> Result<(),FFTError> 
{

    let n = a.len();
    assert_eq!(b.len(), n);

    // first butterfly and twiddle
    b.iter_mut().zip(butterflies.iter()).for_each(|(e, &butterfly)| { *e = a[butterfly]; });
        
    for (butterflies, twiddles) in butterflies[n..].chunks_exact(n)
                                                   .zip(twiddles[n..].chunks_exact(n)) {
    
        // small Fourier transform
        ft_inplace_tf_pow2(b, a);
        
        // butterfly and twiddle
        b.chunks_exact_mut(2).zip(butterflies.chunks_exact(2)
                                             .zip(twiddles.chunks_exact(2)))
                             .for_each(
            |(e, (butterfly, twiddle))| { 
                unsafe {
                    e[0] = *a.get_unchecked(butterfly[0]); 
                    e[1] = *a.get_unchecked(butterfly[1]) * twiddle[1]; 
                }
            });
    }
        
    Ok(()) 
}


pub fn ft_inplace_tf_pow2 (a: &[C], b: &mut [C]) {
    a.chunks_exact(2)
     .zip(b.chunks_exact_mut(2))
     .for_each(|(a_c, b_c)| {
        b_c[0] = a_c[0] + a_c[1];
        b_c[1] = a_c[0] - a_c[1];
    });
}


fn ft_to_ift(y: &mut [C]) {
        
    // rescaling
    let n_r = y.len() as R;
    for a in y.iter_mut() {
        *a /= n_r;
    }

    // swaps
    let n = y.len();
    for i in 1..(n+1)/2 {
        y.swap(i, n-i);
    }
}


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn transpose_1() {
        
        let ft_shape: Vec<usize> = vec![1, 3, 4];
        let ft_struct = MFtStruct::new(ft_shape);

        let x = vec![
            C::new(1.,1.), C::new(2.,1.), C::new(3.,1.), C::new(4.,1.), 
            C::new(1.,2.), C::new(2.,2.), C::new(3.,2.), C::new(4.,2.), 
            C::new(1.,3.), C::new(2.,3.), C::new(3.,3.), C::new(4.,3.), 
        ];

        let mut y = vec![C::new(0.,0.); x.len()];
        ft_struct.transpose(&x, &mut y, 1, 2, 3, 4).unwrap();

        let expected_y = vec![
            C::new(1.,1.), C::new(1.,2.), C::new(1.,3.), 
            C::new(2.,1.), C::new(2.,2.), C::new(2.,3.), 
            C::new(3.,1.), C::new(3.,2.), C::new(3.,3.), 
            C::new(4.,1.), C::new(4.,2.), C::new(4.,3.), 
        ];
        
        assert_eq!(y, expected_y);
    }
    
    #[test]
    fn transpose_2() {
        
        let ft_shape: Vec<usize> = vec![1, 3, 4];
        let ft_struct = MFtStruct::new(ft_shape);

        let x = vec![
            C::new(1.,1.), C::new(2.,1.), C::new(3.,1.), C::new(4.,1.), 
            C::new(1.,2.), C::new(2.,2.), C::new(3.,2.), C::new(4.,2.), 
            C::new(1.,3.), C::new(2.,3.), C::new(3.,3.), C::new(4.,3.), 
        ];

        let mut y = vec![C::new(0.,0.); x.len()];
        ft_struct.transpose(&x, &mut y, 0, 1, 1, 3).unwrap();
        
        assert_eq!(y, x);
    }
}
