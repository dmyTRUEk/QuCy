//! Mathematical Vector (3 dimensional)

use std::ops::{Neg, Sub, Add, Div, Mul};
use crate::math::float;



#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3d {
    pub x: float,
    pub y: float,
    pub z: float
}

impl Vec3d {
    pub fn new(x: float, y: float, z: float) -> Self {
        Self { x, y, z }
    }

    pub fn len2(&self) -> float {
        // self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
        self.x*self.x + self.y*self.y + self.z*self.z
    }

    pub fn len(&self) -> float {
        self.len2().sqrt()
    }

    pub fn dot(&self, rhs: Self) -> float {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }
}

impl Neg for Vec3d {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Vec3d {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Add for Vec3d {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Vec3d {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Sub for Vec3d {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Vec3d {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Div<float> for Vec3d {
    type Output = Self;
    fn div(self, rhs: float) -> Self::Output {
        Vec3d {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl Mul<float> for Vec3d {
    type Output = Self;
    fn mul(self, rhs: float) -> Self::Output {
        Vec3d {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Mul<Vec3d> for float {
    type Output = Vec3d;
    fn mul(self, rhs: Vec3d) -> Self::Output {
        Vec3d {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}

