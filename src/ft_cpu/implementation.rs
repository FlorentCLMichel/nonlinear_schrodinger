use crate::{ R, C, FtStruct, MFtStruct, FFTError };
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
        for d in 0..dimension {
            total_length *= shape[d];
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
    pub fn ft(&self, x: &Vec<C>) -> Result<Vec<C>, FFTError> {
        let mut z = x.clone();
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
    pub fn ift(&self, x: &Vec<C>) -> Result<Vec<C>, FFTError> {
        let mut z = x.clone();
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
            for i in 0..x.len() {
                y[i] = x[i];
            }
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
            for i in 0..x.len() {
                y[i] = x[i];
            }
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
        if i > j {
            let buf = i;
            i = j;
            j = buf;
        }

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
