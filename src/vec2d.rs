//! Mathematical Vector (2 dimensional)

use std::ops::{Neg, Sub, Add, Div, Mul};



#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec2d {
    pub x: f32,
    pub y: f32
}

impl Vec2d {
    pub fn new(x: f32, y: f32) -> Self {
        Vec2d { x, y }
    }

    pub fn len2(&self) -> f32 {
        // self.x.powi(2) + self.y.powi(2)
        self.x*self.x + self.y*self.y
    }

    pub fn len(&self) -> f32 {
        self.len2().sqrt()
    }

    pub fn dot(&self, rhs: Self) -> f32 {
        self.x * rhs.x + self.y * rhs.y
    }

    pub fn is_inside_circle(&self, radius2: f32) -> bool {
        self.len2() < radius2
    }

    pub fn is_inside_rectangle(&self, vec_min: Self, vec_max: Self) -> bool {
        (vec_min.x <= self.x && self.x <= vec_max.x) && (vec_min.y <= self.y && self.y <= vec_max.y)
    }

    pub fn is_inside_triangle(&self, point1: Self, point2: Self, point3: Self) -> bool {
        fn triangle_sign(p1: Vec2d, p2: Vec2d, p3: Vec2d) -> f32 {
            (p1.x-p3.x)*(p2.y-p3.y) - (p2.x-p3.x)*(p1.y-p3.y)
        }
        let d1 = triangle_sign(*self, point1, point2);
        let d2 = triangle_sign(*self, point2, point3);
        let d3 = triangle_sign(*self, point3, point1);
        !((d1 < 0.0 || d2 < 0.0 || d3 < 0.0) && (d1 > 0.0 || d2 > 0.0 || d3 > 0.0))
    }

}

impl Neg for Vec2d {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Vec2d {
            x: -self.x,
            y: -self.y,
        }
    }
}

impl Add for Vec2d {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Vec2d {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Sub for Vec2d {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Vec2d {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl Div<f32> for Vec2d {
    type Output = Self;
    fn div(self, rhs: f32) -> Self::Output {
        Vec2d {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

impl Mul<f32> for Vec2d {
    type Output = Self;
    fn mul(self, rhs: f32) -> Self::Output {
        Vec2d {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Mul<Vec2d> for f32 {
    type Output = Vec2d;
    fn mul(self, rhs: Vec2d) -> Self::Output {
        Vec2d {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}

