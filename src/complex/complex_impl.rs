use crate::complex::*;

impl<T> Complex<T> 
where T: RealNumber
{
    /// Create a complex number from two real ones, interpreted as its real and imaginary parts
    ///
    /// # Arguments
    ///
    /// * `r`: real part
    /// * `i`: imaginary part
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::complex::Complex;
    ///
    /// # fn main() {
    /// #
    /// let r: f64 = 2.;
    /// let i: f64 = -1.;
    ///
    /// let c = Complex::new(r, i);
    ///
    /// assert_eq!(c.real, r);
    /// assert_eq!(c.imag, i);
    /// #
    /// # }
    /// ```
    pub fn new(r: T, i: T) -> Self {
        Complex::<T> {
            real: r,
            imag: i
        }
    }
    
    /// complex conjugation
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::complex::Complex;
    ///
    /// # fn main() {
    /// #
    /// let r: f64 = 2.;
    /// let i: f64 = -1.;
    ///
    /// let c = Complex::new(r, i);
    /// let c_c = c.conjugate();
    ///
    /// assert_eq!(c_c.real, r);
    /// assert_eq!(c_c.imag, -i);
    /// #
    /// # }
    /// ```
    pub fn conjugate(&self) -> Self {
        Complex::<T> {
            real: self.real,
            imag: -self.imag
        }
    }

    /// absolute value
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::complex::Complex;
    ///
    /// # fn main() {
    /// #
    /// let r: f64 = -3.;
    /// let i: f64 = 4.;
    ///
    /// let c = Complex::new(r, i);
    /// let c_abs = c.abs();
    ///
    /// assert_eq!(c_abs, 5.);
    /// #
    /// # }
    /// ```
    pub fn abs(&self) -> T {
        (self.real * self.real + self.imag * self.imag).sqrt()
    }
    
    /// squared absolute value
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::complex::Complex;
    ///
    /// # fn main() {
    /// #
    /// let r: f64 = -3.;
    /// let i: f64 = 4.;
    ///
    /// let c = Complex::new(r, i);
    /// let c_abs2 = c.abs2();
    ///
    /// assert_eq!(c_abs2, 25.);
    /// #
    /// # }
    /// ```
    pub fn abs2(&self) -> T {
        self.real * self.real + self.imag * self.imag
    }
    
    /// complex exponentiation
    ///
    /// # Example
    ///
    /// ```
    /// use nonlinear_schrodinger::complex::Complex;
    ///
    /// # fn main() {
    /// #
    /// let r: f64 = 0.;
    /// let i: f64 = 0.;
    ///
    /// let c = Complex::new(r, i);
    /// let c_exp = c.exp();
    ///
    /// assert_eq!(c_exp, Complex::new(1.,0.));
    /// #
    /// # }
    /// ```
    pub fn exp(&self) -> Self {
        Complex::<T> {
            real: self.imag.cos(),
            imag: self.imag.sin()
        } * self.real.exp()
    }
}


impl<T: RealNumber> std::cmp::PartialEq for Complex<T> {
    fn eq(&self, other: &Self) -> bool {
        (self.real) == (other.real) && (self.imag) == (other.imag)
    }
}

impl<T: RealNumber> ops::Add for Complex<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Complex::<T> {
            real: self.real + rhs.real,
            imag: self.imag + rhs.imag
        }
    }
}

impl<T: RealNumber> ops::AddAssign for Complex<T> {
    fn add_assign(&mut self, rhs: Self) {
        self.real += rhs.real;
        self.imag += rhs.imag;
    }
}

impl<T: RealNumber> ops::Sub for Complex<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Complex::<T> {
            real: self.real - rhs.real,
            imag: self.imag - rhs.imag
        }
    }
}

impl<T: RealNumber> ops::SubAssign for Complex<T> {
    fn sub_assign(&mut self, rhs: Self) {
        self.real -= rhs.real;
        self.imag -= rhs.imag;
    }
}

impl<T: RealNumber> ops::Mul<Self> for Complex<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Complex::<T> {
            real: self.real * rhs.real - self.imag * rhs.imag,
            imag: self.real * rhs.imag + self.imag * rhs.real
        }
    }
}

impl<T: RealNumber> ops::MulAssign<Self> for Complex<T> {
    fn mul_assign(&mut self, rhs: Self) {
        let r = self.real;
        self.real = self.real * rhs.real - self.imag * rhs.imag;
        self.imag = r * rhs.imag + self.imag * rhs.real;
    }
}

impl<T: RealNumber> ops::Mul<T> for Complex<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        Complex::<T> {
            real: self.real * rhs,
            imag: self.imag * rhs
        }
    }
}

impl<T: RealNumber> ops::MulAssign<T> for Complex<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.real *= rhs;
        self.imag *= rhs;
    }
}

impl<T: RealNumber> ops::Div<T> for Complex<T> {
    type Output = Self;
    fn div(self, rhs: T) -> Self {
        Complex::<T> {
            real: self.real / rhs,
            imag: self.imag / rhs
        }
    }
}

impl<T: RealNumber> ops::Div<Self> for Complex<T> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self * rhs.conjugate() / (self.abs() * rhs.abs())
    }
}

impl<T: RealNumber> ops::DivAssign<T> for Complex<T> {
    fn div_assign(&mut self, rhs: T) {
        self.real /= rhs;
        self.imag /= rhs;
    }
}

impl<T: RealNumber> ops::DivAssign<Self> for Complex<T> {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self * rhs.conjugate() / ((*self).abs() * rhs.abs())
    }
}

impl<T: RealNumber> std::fmt::Display for Complex<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.imag >= T::zero() {
            write!(f, "{}+{}i", self.real, self.imag)
        } else {
            write!(f, "{}{}i", self.real, self.imag)
        }
    }
}

impl<T: RealNumber> std::fmt::Debug for Complex<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Complex")
         .field("real", &self.real)
         .field("imag", &self.imag)
         .finish()
    }
}
