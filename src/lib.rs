extern crate nalgebra as na;
use core::cell::{BorrowError, BorrowMutError, Ref, RefCell};
use na::{point, vector, Point3, Vector3};

pub const G: f64 = 0.0376915586937436; // Units: Jupiter Mass - AU - Years
pub const SOFTENER: f64 = 0.01; // Artificially decrease force strength relative to seperation

fn specific_acceleration(seperation: Vector3<f64>) -> Vector3<f64> {
    (G / (seperation.magnitude_squared() + SOFTENER)) * seperation.normalize()
}

pub fn test_calc(test: &TestParticle, massive: &MassiveParticle) -> Result<(), BorrowMutError> {
    test.update_acceleration(massive, &specific_acceleration(test.seperation(massive)))
}

pub fn massive_calc(a: &MassiveParticle, b: &MassiveParticle) -> Result<(), BorrowMutError> {
    let specific_acceleration = &specific_acceleration(a.seperation(b));
    a.update_acceleration(b, &specific_acceleration)?;
    b.update_acceleration(a, &-specific_acceleration)
}

pub fn particle_compute(a: &impl Particle, b: &impl Particle) -> Result<(), BorrowMutError> {
    let specific_acceleration = &specific_acceleration(a.seperation(b));
    a.update_acceleration(b, &specific_acceleration)?;
    b.update_acceleration(a, &-specific_acceleration)
}

pub trait Particle {
    fn mass(&self) -> f64;
    fn position(&self) -> &Point3<f64>;
    fn velocity(&self) -> &Vector3<f64>;
    fn acceleration(&self) -> Result<Ref<Vector3<f64>>, BorrowError>;
    fn update(&mut self, integration_interval: f64) -> Result<(), BorrowError>;

    fn seperation(&self, other: &impl Particle) -> Vector3<f64> {
        other.position() - self.position()
    }

    fn store_acceleration(&self, acceleration: &Vector3<f64>) -> Result<(), BorrowMutError>;

    fn update_acceleration(
        &self,
        other: &impl Particle,
        specific_acceleration: &Vector3<f64>,
    ) -> Result<(), BorrowMutError> {
        self.store_acceleration(&(specific_acceleration * other.mass()))
    }

    fn new_position(&self, integration_interval: f64) -> Result<Point3<f64>, BorrowError> {
        Ok(self.position()
            + (*(self.acceleration()?) * (integration_interval.powi(2) / 2.0)
                + self.velocity() * integration_interval))
    }

    fn new_velocity(&self, integration_interval: f64) -> Result<Vector3<f64>, BorrowError> {
        Ok((*(self.acceleration()?) * integration_interval) + self.velocity())
    }
}

#[derive(PartialEq, Clone, Default, Debug)]
pub struct TestParticle {
    position: Point3<f64>,
    velocity: Vector3<f64>,
    acceleration: RefCell<Vector3<f64>>,
}

impl Particle for TestParticle {
    fn mass(&self) -> f64 {
        0.0
    }

    fn position(&self) -> &Point3<f64> {
        &self.position
    }

    fn velocity(&self) -> &Vector3<f64> {
        &self.velocity
    }

    fn acceleration(&self) -> Result<Ref<Vector3<f64>>, BorrowError> {
        self.acceleration.try_borrow()
    }

    fn update(&mut self, integration_interval: f64) -> Result<(), BorrowError> {
        self.position = self.new_position(integration_interval)?;
        self.velocity = self.new_velocity(integration_interval)?;
        Ok(())
    }

    fn store_acceleration(&self, acceleration: &Vector3<f64>) -> Result<(), BorrowMutError> {
        *self.acceleration.try_borrow_mut()? += acceleration;
        Ok(())
    }
}

impl TestParticle {
    pub fn new(position: Point3<f64>, velocity: Vector3<f64>) -> Self {
        Self {
            position,
            velocity,
            acceleration: RefCell::new(Vector3::<f64>::default()),
        }
    }

    pub const fn new_const(x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) -> Self {
        Self {
            position: point![x, y, z],
            velocity: vector![vx, vy, vz],
            acceleration: RefCell::new(vector![0.0, 0.0, 0.0]),
        }
    }
}

#[derive(PartialEq, Clone, Default, Debug)]
pub struct MassiveParticle {
    mass: f64,
    centre: TestParticle,
}

impl Particle for MassiveParticle {
    fn mass(&self) -> f64 {
        self.mass
    }

    fn position(&self) -> &Point3<f64> {
        self.centre.position()
    }

    fn velocity(&self) -> &Vector3<f64> {
        self.centre.velocity()
    }

    fn acceleration(&self) -> Result<Ref<Vector3<f64>>, BorrowError> {
        self.centre.acceleration()
    }

    fn update(&mut self, integration_interval: f64) -> Result<(), BorrowError> {
        self.centre.update(integration_interval)
    }

    fn store_acceleration(&self, acceleration: &Vector3<f64>) -> Result<(), BorrowMutError> {
        self.centre.store_acceleration(acceleration)
    }
}

impl MassiveParticle {
    pub fn new(position: Point3<f64>, velocity: Vector3<f64>, mass: f64) -> Self {
        Self {
            mass,
            centre: TestParticle {
                position,
                velocity,
                acceleration: RefCell::new(Vector3::<f64>::default()),
            },
        }
    }

    pub fn from_test(centre: TestParticle, mass: f64) -> Self {
        Self { mass, centre }
    }

    pub const fn new_const(x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64, m: f64) -> Self {
        Self {
            mass: m,
            centre: TestParticle::new_const(x, y, z, vx, vy, vz),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    use itertools::Itertools;

    const P0V0M0: TestParticle = TestParticle::new_const(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    const P1V0M0: TestParticle = TestParticle::new_const(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    const P0V0M1: MassiveParticle = MassiveParticle::new_const(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

    const P1V0M1: MassiveParticle = MassiveParticle::new_const(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

    #[test]
    fn seperation() {
        let seperation = P0V0M0.seperation(&P1V0M0);
        assert_eq!(seperation, vector![1.0, 0.0, 0.0]);
    }

    #[test]
    fn specific_acceleration() {
        fn specific_acceleration_test(seperation: &Vector3<f64>)  -> Vector3<f64> {
            // (G / (seperation.magnitude_squared() + SOFTENER)) * seperation.normalize()
            //mut (G / (seperation.normalize_mut().powi(2) + SOFTENER)) * seperation
            G * (*seperation / (seperation.magnitude_squared().powf(1.5) + SOFTENER))
        }

        for nums in (-1000000..1000000).step_by(10000).combinations(3) {
            let v = vector![f64::from(nums[0]), f64::from(nums[1]), f64::from(nums[2])];

            let t = specific_acceleration_test(&v);
            let f = crate::specific_acceleration(v);
            let d = f.cross(&t);

            assert!(d.magnitude() <= 0.00000000000000000001, "\nMagnitude: {:?}\nVector tested: {:?}\nlib:  {:?}\ntest: {:?}\n", d.magnitude(), v, f, t);
        }
    }

    #[test]
    fn one_iteration() {
        let dt = 0.0001;
        let dT = 4.0;
        let T = 12.0;

        const VMAX: f64 = 0.453;
        let mut V = VMAX;

        let mut v1 = crate::MassiveParticle::new_const(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1047.57);
        let mut v2 = crate::MassiveParticle::new_const(5.035, 0.0, 0.0, 0.0, V, 0.0, 1.0);

        println!("{:#?}\n{:#?}", v1, v2);

        for t in (0..dT as u64) {
            for t in (0..(T/(dt*dT)) as u64) {
                crate::particle_compute(&v1, &v2);
                v1.update(dt);
                v2.update(dt);

                V = V.max(v2.velocity().magnitude());
            }

            println!("{:#?}\n{:#?}", v1, v2);
        }

            assert!(V < VMAX*1.1);
    }

    // #[test]
    // fn acceleration() {

    // }

    // #[test]
    // fn update_acceleration() {
    //     let t = p1v0m0.clone();

    //     t.update_acceleration(p0v0m1);
    // }
}
