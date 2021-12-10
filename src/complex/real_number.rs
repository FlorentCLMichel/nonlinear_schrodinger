use std::ops;

/// A trait that must be satisfied by ‘acceptable’ real number types
pub trait RealNumber : Copy + Clone + ops::Add<Output=Self> 
                                    + ops::Sub<Output=Self> 
                                    + ops::Mul<Output=Self> 
                                    + ops::Div<Output=Self> 
                                    + ops::Neg<Output=Self> 
                                    + ops::AddAssign 
                                    + ops::SubAssign 
                                    + ops::MulAssign 
                                    + ops::DivAssign 
                                    + std::cmp::PartialEq
                                    + std::cmp::PartialOrd
                                    + std::fmt::Display
                                    + std::fmt::Debug
                                    + std::convert::From<i8>
{
    fn exp(self) -> Self;
    fn cos(self) -> Self;
    fn sin(self) -> Self;
    fn sqrt(self) -> Self;
    fn zero() -> Self;
}

#[macro_export]
macro_rules! arithmetic_impl {
    ($T: ident) => (
    impl RealNumber for $T
    {
        fn exp(self) -> $T {
            $T::exp(self)
        }
        fn cos(self) -> $T {
            $T::cos(self)
        }
        fn sin(self) -> $T {
            $T::sin(self)
        }
        fn sqrt(self) -> $T {
            $T::sqrt(self)
        }
        fn zero() -> $T { 0. }
    }
)}


// implement the RealNumber trait for f32 and f64
arithmetic_impl!(f32);
arithmetic_impl!(f64);
