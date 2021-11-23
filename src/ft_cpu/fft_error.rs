//! Defines a generic error for Fourier transforms and related computations

#[derive(Debug)]
/// Error for the fft and related functions
pub struct FFTError {
    message: Box<String>,
}


impl std::fmt::Display for FFTError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "FFTError: {}", *self.message)
    }
}


impl std::error::Error for FFTError {}


impl FFTError {

    /// create a new `FFTError`
    ///
    /// # Argument
    ///
    /// * `s`: the error message
    pub fn new(s: String) -> FFTError {
        FFTError {message: Box::new(s)}
    }
}
