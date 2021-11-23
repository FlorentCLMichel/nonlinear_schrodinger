use crate::{ R, C, FtStruct, FFTError };
use crate::ft_cpu::plan::*;

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

        // find the optimal plan
        let plan = find_plan(n);

        // find the butterfly coefficients and twiddle factors
        let (butterflies, twiddles, twiddles_small_ft) = get_butterflies_and_twiddles(n, &plan);

        // return the ft structure
        FtStruct { n, plan, butterflies, twiddles, twiddles_small_ft }
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

        let mut first_index_small_ft = 0;
        for k in 0..self.plan.len() {
    
            // butterfly and twiddle
            for i in 0..self.n {
                b[i] = a[self.butterflies[k*self.n+i]] * self.twiddles[k*self.n+i];
            }
        
            // small Fourier transform
            ft_inplace_tf(b, a, self.plan[self.plan.len()-1-k], 
                          &self.twiddles_small_ft[first_index_small_ft..
                                        first_index_small_ft + self.plan[self.plan.len()-1-k]]);
            first_index_small_ft += self.plan[self.plan.len()-1-k];
    
        }
            
        // final butterfly and twiddle
        let k = self.plan.len();
        for i in 0..self.n {
            b[i] = a[self.butterflies[k*self.n+i]] * self.twiddles[k*self.n+i];
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
