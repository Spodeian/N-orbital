extern crate nalgebra as na;
use na::{Point3, Vector3};
// use serde::{Serialize, Deserialize};

// use physical_constants::NEWTONIAN_CONSTANT_OF_GRAVITATION;

// pub const G: f64 = NEWTONIAN_CONSTANT_OF_GRAVITATION;

pub const G: f64 = 0.0376915586937436; // Units: Jupiter Mass - AU - Years
pub const SOFTENER: f64 = 0.0001; // Artificially decrease force strength relative to seperation

pub trait Particle {
    fn mass(&self) -> f64;
    fn position(&self) -> Point3<f64>;
    fn velocity(&self) -> Vector3<f64>;
    fn update(&mut self, acceleration: &Vector3<f64>, integration_interval: f64);

    // new methods
    fn seperation(&self, other: &impl Particle) -> Vector3<f64> {
        other.position() - self.position()
    }

    fn specific_acceleration(&self, mut seperation: Vector3<f64>) -> Vector3<f64> {
        (G / (seperation.normalize_mut().powi(2) + SOFTENER)) * seperation
    }

    fn store_acceleration(&mut self, acceleration: &Vector3<f64>);

    fn update_acceleration(&mut self, other: &impl Particle, specific_acceleration: &Vector3<f64>) {
        self.store_acceleration(&(specific_acceleration * other.mass()))
    }

    // old methods
    fn accel_towards(&self, other: &impl Particle) -> Vector3<f64> {
        let seperation: Vector3<f64> = other.position() - self.position();
        (G * other.mass() / (seperation.magnitude_squared() + SOFTENER)) * (seperation.normalize())
    }

    fn new_position(&self, acceleration: &Vector3<f64>, integration_interval: f64) -> Point3<f64> {
        self.position()
            + (acceleration * (integration_interval.powi(2) / 2.0)
                + self.velocity() * integration_interval)
    }

    fn new_velocity(&self, acceleration: &Vector3<f64>, integration_interval: f64) -> Vector3<f64> {
        (acceleration * integration_interval) + self.velocity()
    }
}

#[derive(PartialEq, Clone, Default, Debug)]
pub struct TestParticle {
    position: Point3<f64>,
    velocity: Vector3<f64>,
    acceleration: Vector3<f64>,
}

impl Particle for TestParticle {
    fn mass(&self) -> f64 {
        0.0
    }

    fn position(&self) -> Point3<f64> {
        self.position
    }

    fn velocity(&self) -> Vector3<f64> {
        self.velocity
    }

    fn update(&mut self, acceleration: &Vector3<f64>, integration_interval: f64) {
        self.position = self.new_position(acceleration, integration_interval);
        self.velocity = self.new_velocity(acceleration, integration_interval);
    }

    fn store_acceleration(&mut self, acceleration: &Vector3<f64>) {
        self.acceleration += acceleration;
    }
}

impl TestParticle {
    pub fn new(position: Point3<f64>, velocity: Vector3<f64>) -> Self {
        Self {
            position,
            velocity,
            acceleration: Vector3::<f64>::default(),
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

    fn position(&self) -> Point3<f64> {
        self.centre.position
    }

    fn velocity(&self) -> Vector3<f64> {
        self.centre.velocity
    }

    fn update(&mut self, acceleration: &Vector3<f64>, integration_interval: f64) {
        self.centre.position = self.new_position(acceleration, integration_interval);
        self.centre.velocity = self.new_velocity(acceleration, integration_interval);
    }

    fn store_acceleration(&mut self, acceleration: &Vector3<f64>) {
        self.centre.store_acceleration(acceleration);
    }
}

impl MassiveParticle {
    pub fn new(position: Point3<f64>, velocity: Vector3<f64>, mass: f64) -> Self {
        Self {
            mass,
            centre: TestParticle {
                position,
                velocity,
                acceleration: Vector3::<f64>::default(),
            },
        }
    }

    pub fn from_test(centre: TestParticle, mass: f64) -> Self {
        Self { mass, centre }
    }
}
