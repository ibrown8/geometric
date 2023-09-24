use num_traits::{Num, float::Float};
use core::ops::{Add, Sub, Mul, Div, AddAssign, MulAssign, DivAssign, SubAssign};
use core::convert::From;
use core::fmt::{Debug, Formatter, Result};
//Dot product
pub trait Dot<RHS> {
    type Output;
    fn dot(&self, other : &RHS) -> Self::Output;
}
//Wedge product
pub trait Wedge<RHS> {
    type Output;
    fn wedge(&self, other : &RHS) -> Self::Output;
}
//Geometric Product
pub trait Geometric<RHS> {
    type Output;
    fn geometric(&self, other : &RHS) -> Self::Output;
}
pub trait Cross {
    type Output;
    fn cross(&self, other : &Self) -> Self::Output;
}
#[derive(Copy, Clone)]
pub struct Vec2<T: Copy> {
    pub x : T,
    pub y : T,
}

macro_rules! impl_bin_op_vector2 {
    ($t:ident, $f:ident) => {
        impl<T: Num + Copy> $t<Vec2<T>> for Vec2<T> {
            type Output = Self;
            #[inline]
            fn $f(self, other : Self) -> Self::Output {
                Vec2 {
                    x : <T as $t>::$f(self.x, other.x),
                    y : <T as $t>::$f(self.y, other.y)
                }
            }
        }
        impl<'a, T: Num + Copy> $t<&'a Vec2<T>> for &'a Vec2<T> {
            type Output = Vec2<T>;
            #[inline]
            fn $f(self, other : Self) -> Self::Output {
                Vec2 {
                    x : <T as $t>::$f(self.x, other.x),
                    y : <T as $t>::$f(self.y, other.y)
                }
            }
        }
        impl<T: Num + Copy> $t<T> for Vec2<T> {
            type Output = Self;
            #[inline]
            fn $f(self, other : T) -> Self::Output {
                Vec2 {
                    x : <T as $t>::$f(self.x, other),
                    y : <T as $t>::$f(self.y, other)
                }
            }
        }
        impl<'a, T: Num + Copy> $t<T> for &'a Vec2<T> {
            type Output = Vec2<T>;
            #[inline]
            fn $f(self, other : T) -> Self::Output {
                Vec2 {
                    x : <T as $t>::$f(self.x, other),
                    y : <T as $t>::$f(self.y, other)
                }
            }
        }
    }
}

macro_rules! impl_bin_op_assign_vector2 {
    ($t:ident, $f:ident) => {
        impl<T: $t + Num + Copy> $t<Vec2<T>> for Vec2<T> {
            #[inline]
            fn $f(&mut self, other : Self){
                <T as $t>::$f(&mut self.x, other.x);
                <T as $t>::$f(&mut self.y, other.y);
            }
        }
        impl<T: $t + Num + Copy> $t<T> for Vec2<T> {
            #[inline]
            fn $f(&mut self, other : T){
                <T as $t>::$f(&mut self.x, other);
                <T as $t>::$f(&mut self.y, other);
            }
        }
    }
}
#[inline]
pub fn lerp<X, T>(v0 : X, v1 : X, t : T) -> X 
    where 
        T: Float + Copy,
        X: Add<Output = X> + Mul<T, Output = X>
{
      v0 * (T::one() - t) + v1 * t
}
impl_bin_op_vector2!{Add, add}
impl_bin_op_vector2!{Sub, sub}
impl_bin_op_vector2!{Mul, mul}
impl_bin_op_vector2!{Div, div}
impl_bin_op_assign_vector2!{AddAssign, add_assign}
impl_bin_op_assign_vector2!{SubAssign, sub_assign}
impl_bin_op_assign_vector2!{MulAssign, mul_assign}
impl_bin_op_assign_vector2!{DivAssign, div_assign}

impl <T: Num + Copy> Dot<Vec2<T>> for Vec2<T> {
    type Output = T;
    #[inline]
    fn dot(&self, other : &Self) -> T {
        self.x * other.x + self.y * other.y
    }
}
//Wedge product
impl <T: Num + Copy> Wedge<Vec2<T>> for Vec2<T> {
    type Output = T;
    #[inline]
    fn wedge(&self, other : &Self) -> T {
        self.x * other.y - self.y * other.x
    }
}

impl<T: Num + Copy + Debug> Debug for Vec2<T> {
    fn fmt(&self, f : &mut Formatter<'_>) -> Result {
        f.debug_struct("Vec2")
            .field("x", &self.x)
            .field("y", &self.y)
            .finish()
    }
}
 
pub struct Rotor2<T: Copy + Num> {
    pub scalar : T,
    pub pseudo : T
}

impl <T: Num + Copy> Geometric<Vec2<T>> for Vec2<T> {
    type Output = Rotor2<T>;
    #[inline]
    fn geometric(&self, other : &Self) -> Self::Output {
        Rotor2 {
            scalar : self.dot(other),
            pseudo : self.wedge(other)
        }
    }
}

impl <T: Num + Copy> Geometric<Rotor2<T>> for Vec2<T> {
    type Output = Vec2<T>;
    #[inline]
    fn geometric(&self, other : &Rotor2<T>) -> Self::Output {
        Vec2 {
            x : self.x * other.scalar - self.y * other.pseudo,
            y : self.x * other.pseudo + self.y * other.scalar
        }
    }
}

impl <T: Num + Copy> Geometric<Vec2<T>> for Rotor2<T> {
    type Output = Vec2<T>;
    #[inline]
    fn geometric(&self, other : &Vec2<T>) -> Self::Output {
        Vec2 {
            x : self.scalar * other.x + self.pseudo * other.y,
            y : self.scalar * other.y - self.pseudo * other.x
        }
    }
}

impl <T: Num + Copy> Geometric<Rotor2<T>> for Rotor2<T> {
    type Output = Rotor2<T>;
    #[inline]
    fn geometric(&self, other : &Self) -> Self::Output { 
        Rotor2 {
            scalar : self.scalar * other.scalar - self.pseudo * other.pseudo,
            pseudo : self.scalar * other.pseudo + self.pseudo * other.scalar
        }
    }
}

impl<T: Float> Rotor2<T> {
    //Euler identity for Bivectors
    #[inline]
    pub fn from_angle(theta : T) -> Self {
        let (scalar, pseudo) = theta.sin_cos();
        Self {scalar, pseudo}
    }
}
impl<T: Float> Vec2<T> {
    #[inline]
    pub fn len_squared(&self) -> T {
        self.dot(self)
    }
    #[inline]
    pub fn len(&self) -> T {
        self.x.hypot(self.y)
    }
    #[inline]
    pub fn normalize(&self) -> Self {
        let inv_len = T::one() / self.len();
        self * inv_len
    }
}
#[derive(Copy, Clone)]
pub struct Vec3<T: Copy> {
    pub x : T,
    pub y : T,
    pub z : T,
}

macro_rules! impl_bin_op_vector3 {
    ($t:ident, $f:ident) => {
        impl<T: Num + Copy> $t<Vec3<T>> for Vec3<T> {
            type Output = Vec3<T>;
            #[inline]
            fn $f(self, other : Self) -> Self::Output {
                Vec3 {
                    x : <T as $t>::$f(self.x, other.x),
                    y : <T as $t>::$f(self.y, other.y),
                    z : <T as $t>::$f(self.z, other.z)
                }
            }
        }
        impl<'a, T: Num + Copy> $t<&'a Vec3<T>> for &'a Vec3<T> {
            type Output = Vec3<T>;
            #[inline]
            fn $f(self, other : Self) -> Self::Output {
                Vec3 {
                    x : <T as $t>::$f(self.x, other.x),
                    y : <T as $t>::$f(self.y, other.y),
                    z : <T as $t>::$f(self.z, other.z)
                }
            }
        }
        impl<T: Num + Copy> $t<T> for Vec3<T> {
            type Output = Self;
            #[inline]
            fn $f(self, other : T) -> Self::Output {
                Vec3 {
                    x : <T as $t>::$f(self.x, other),
                    y : <T as $t>::$f(self.y, other),
                    z : <T as $t>::$f(self.z, other)
                }
            }
        }
        impl<'a, T: Num + Copy> $t<T> for &'a Vec3<T> {
            type Output = Vec3<T>;
            #[inline]
            fn $f(self, other : T) -> Self::Output {
                Vec3 {
                    x : <T as $t>::$f(self.x, other),
                    y : <T as $t>::$f(self.y, other),
                    z : <T as $t>::$f(self.z, other)
                }
            }
        }
    }
}

macro_rules! impl_bin_op_assign_vector3 {
    ($t:ident, $f:ident) => {
        impl<T: $t + Num + Copy> $t<Vec3<T>> for Vec3<T> {
            #[inline]
            fn $f(&mut self, other : Self){
                <T as $t>::$f(&mut self.x, other.x);
                <T as $t>::$f(&mut self.y, other.y);
                <T as $t>::$f(&mut self.z, other.z);
            }
        }
        impl<T: $t + Num + Copy> $t<T> for Vec3<T> {
            #[inline]
            fn $f(&mut self, other : T){
                <T as $t>::$f(&mut self.x, other);
                <T as $t>::$f(&mut self.y, other);
                <T as $t>::$f(&mut self.z, other);
            }
        }
    }
}

impl_bin_op_vector3!{Add, add}
impl_bin_op_vector3!{Sub, sub}
impl_bin_op_vector3!{Mul, mul}
impl_bin_op_vector3!{Div, div}
impl_bin_op_assign_vector3!{AddAssign, add_assign}
impl_bin_op_assign_vector3!{SubAssign, sub_assign}
impl_bin_op_assign_vector3!{MulAssign, mul_assign}
impl_bin_op_assign_vector3!{DivAssign, div_assign}

impl<T: Num + Copy + Debug> Debug for Vec3<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        f.debug_struct("Vec3")
            .field("x", &self.x)
            .field("y", &self.y)
            .field("z", &self.z)
            .finish()
    }
}
//Wedge product
impl <T: Num + Copy> Dot<Vec3<T>> for Vec3<T> {
    type Output = T;
    #[inline]
    fn dot(&self, other : &Self) -> Self::Output {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}

impl<T: Float> Vec3<T> {
    #[inline]
    pub fn len_squared(&self) -> T {
        self.dot(self)
    }
    pub fn len(&self) -> T {
        self.len_squared().sqrt()
    }
    #[inline]
    pub fn normalize(&self) -> Self {
        let inv_len = T::one() / self.len();
        self * inv_len
    }
}
#[derive(Copy, Clone)]
pub struct Bivector3<T: Copy + Num> {
    pub xy : T,
    pub yz : T,
    pub zx : T
}

impl<T: Num + Copy> From<Bivector3<T>> for Vec3<T> {
    fn from(value : Bivector3<T>) -> Self {
        Vec3 {x : value.yz, y : value.zx, z : value.xy}
    } 
}

impl<T: Num + Copy> Add<Bivector3<T>> for Bivector3<T> {
    type Output = Bivector3<T>;
    #[inline]
    fn add(self, other : Self) -> Self::Output {
        Bivector3 {
            xy : self.xy + other.xy,
            yz : self.yz + other.yz,
            zx : self.zx + other.zx
        }
    }
}

impl<T: Num + Copy> Sub<Bivector3<T>> for Bivector3<T> {
    type Output = Bivector3<T>;
    #[inline]
    fn sub(self, other : Self) -> Self::Output{
        Bivector3 {
            xy : self.xy - other.xy,
            yz : self.yz - other.yz,
            zx : self.zx - other.zx
        }
    }
}

impl<T: Num + Copy> Mul<T> for  Bivector3<T> {
    type Output = Bivector3<T>;
    #[inline]
    fn mul(self, other : T) -> Self::Output {
        Bivector3 {
            xy : self.xy * other,
            yz : self.yz * other,
            zx : self.zx * other
        }    
    }
}

impl<'a, T: Num + Copy> Mul<T> for &'a Bivector3<T> {
    type Output = Bivector3<T>;
    #[inline]
    fn mul(self, other : T) -> Self::Output {
        Bivector3 {
            xy : self.xy * other,
            yz : self.yz * other,
            zx : self.zx * other
        }    
    }
}

impl<T: Num + Copy> Dot<Bivector3<T>> for Bivector3<T> {
    type Output = T;
    #[inline]
    fn dot(&self, other : &Self) -> Self::Output {
        T::zero() - (self.xy * other.xy) - (self.yz * other.yz) - (self.zx * other.zx)
    }
}
//Wedge product
impl <'a, T: Num + Copy> Wedge<Vec3<T>> for Vec3<T> {
    type Output = Bivector3<T>;
    #[inline]
    fn wedge(&self, other : &Self) -> Self::Output {
        Bivector3{
            xy : self.x * other.y - self.y * other.x,
            yz : self.y * other.z - self.z * other.y,
            zx : self.z * other.x - self.x * other.z
        }
    }
}
impl<T: Copy + Num> Cross for Vec3<T> {   
    type Output = Self;
    #[inline]
    fn cross(&self, other : &Self) -> Self {
        Self::from(self.wedge(other))
    }
}
impl<T: Num + Copy> Dot<Bivector3<T>> for Vec3<T> {
    type Output = Vec3<T>;
    #[inline]
    fn dot(&self, other : &Bivector3<T>) -> Self::Output {
        Vec3 {
            x : self.z * other.zx - self.y * other.xy,
            y : self.x * other.xy - self.z * other.yz,
            z : self.y * other.yz - self.x * other.zx
        }
    }
}
impl<T: Num + Copy> Cross for Bivector3<T> {
    type Output = Self;
    #[inline]
    fn cross(&self, other : &Self) -> Self {
        Self {
            xy : self.yz * other.zx - self.zx * other.yz,
            yz : self.zx * other.xy - self.xy * other.zx,
            zx : self.xy * other.yz - self.yz * other.xy
        }
    }
}
//Quaternion
pub struct Rotor3<T: Copy + Num> {
    pub scalar : T,
    pub bivec : Bivector3<T>
}
impl<T: Num + Copy> Vec3<T> {
    #[inline]
    //Geometric product
    pub fn geometric(&self, other : &Self) -> Rotor3<T> {
        Rotor3 {
            scalar : self.dot(other),
            bivec : self.wedge(other)
        }
    }
}

impl<T: Num + Copy> Geometric<Rotor3<T>> for Rotor3<T> {
    type Output = Rotor3<T>;
    fn geometric(&self, other : &Self) -> Self::Output {
        Self {
            scalar : self.scalar * other.scalar + self.bivec.dot(&other.bivec),
            bivec : (&other.bivec * self.scalar) + (&self.bivec * other.scalar) + self.bivec.cross(&other.bivec)
        }
    }
}
//TODO: Fixme!
//https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
impl Rotor3<f32> {
    pub fn from_euler(roll : f32, pitch : f32, yaw : f32) -> Self {
        let (sr, cr) = (roll / 2.0).sin_cos();
        let (sp, cp) = (pitch / 2.0).sin_cos();
        let (sy, cy) = (yaw / 2.0).sin_cos();
        let w = (cr * cp * cy) + (sr * sp * sy);
        let x = (sr * cp * cy) - (cr * sp * sy);
        let y = (cr * sp * cy) + (sr * cp * sy);
        let z = (cr * cp * sy) - (sr * sp * cy);
        Rotor3 {
            scalar : w,
            bivec : Bivector3 {
                xy : z, yz : x, zx : y
            }
        }
    }
   /* pub fn to_matrix4(&self) -> Matrix4<f32> {
        let w = self.scalar;
        let x = self.bivec.yz;
        let y = self.bivec.zx;
        let z = self.bivec.xy;
        let m_00 = 1.0 - (2.0 * ((y * y) + (z * z)));
        let m_11 = 1.0 - (2.0 * ((z * z) + (x * x)));
        let m_22 = 1.0 - (2.0 * ((x * x) + (y * y)));
        let m_01 = 2.0 * ((x * y) - (w * z));
        let m_10 = 2.0 * ((x * y) + (w * z));
        let m_20 = 2.0 * ((z * x) - (w * y));
        let m_02 = 2.0 * ((z * x) + (w * y));
        let m_12 = 2.0 * ((y * z) - (w * x));
        let m_21 = 2.0 * ((y * z) + (w * x));
        let result = Matrix4 {
            cols : [[m_00, m_10, m_20, 0.0], 
                    [m_01, m_11, m_21, 0.0], 
                    [m_02, m_12, m_22, 0.0],
                    [0.0, 0.0, 0.0, 1.0]]
        };
        result
    } */
}
//col major
struct Matrix3<T: Num + Copy> {
    pub cols : [Vec3<T>; 3]
}
impl From<Rotor3<f32>> for Matrix3<f32> {
    fn from(rot : Rotor3<f32>) -> Self {
        let w = rot.scalar;
        let x = rot.bivec.yz;
        let y = rot.bivec.zx;
        let z = rot.bivec.xy;
        let m_00 = 1.0 - (2.0 * ((y * y) + (z * z)));
        let m_11 = 1.0 - (2.0 * ((z * z) + (x * x)));
        let m_22 = 1.0 - (2.0 * ((x * x) + (y * y)));
        let m_01 = 2.0 * ((x * y) - (w * z));
        let m_10 = 2.0 * ((x * y) + (w * z));
        let m_20 = 2.0 * ((z * x) - (w * y));
        let m_02 = 2.0 * ((z * x) + (w * y));
        let m_12 = 2.0 * ((y * z) - (w * x));
        let m_21 = 2.0 * ((y * z) + (w * x));
        Self {
            cols : [Vec3{x : m_00, y : m_10, z : m_20},
                    Vec3{x : m_01, y : m_11, z : m_21},
                    Vec3{x : m_02, y : m_12, z : m_22}
            ]
        }
    }
}
impl<'a, T: Num + Copy> Mul<Vec3<T>> for &'a Matrix3<T> {
    type Output = Vec3<T>;
    fn mul(self, other : Vec3<T>) -> Self::Output {
        (self.cols[0] * other.x) + (self.cols[1] * other.y) + (self.cols[2] * other.z)
    }
}
pub struct Vec4<T: Copy + Num> {
    pub x : T,
    pub y : T,
    pub z : T,
    pub w : T
}

pub struct Matrix4<T: Copy + Num> {
    pub cols : [Vec4<T>; 4]
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
