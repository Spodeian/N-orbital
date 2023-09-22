use csv::WriterBuilder;
use nalgebra::Vector3;
use std::error::Error;
use std::time::Instant;
use N_orbital::*;

fn main() -> Result<(), Box<dyn Error>> {
    let start = Instant::now();

    // Setup important timers
    const INTEGRATION_INTERVAL: f64 = 0.001; // Exact
    let mut output_interval: f64 = 0.01; // Approximate
    let mut total_time: f64 = 100.0; // Approximate

    println!("Integration Inverval: {:?}", INTEGRATION_INTERVAL);
    println!("Input Output Inverval: {:?}", output_interval);
    println!("Input Total Time: {:?}", total_time);

    // Normalise
    output_interval = (output_interval / INTEGRATION_INTERVAL).floor() * INTEGRATION_INTERVAL;
    total_time = (total_time / output_interval).floor() * output_interval;

    println!("Calculated Output Inverval: {:?}", output_interval);
    println!("Calculated Total Time: {:?}", total_time);

    // Some test data
    let sun = MassiveParticle::new([0.0; 3].into(), [0.0; 3].into(), 1048.0);
    let mercury =
        MassiveParticle::new([0.4, 0.0, 0.0].into(), [0.0, 10.020, 0.0].into(), 0.0001739);
    let venus = MassiveParticle::new([0.723, 0.0, 0.0].into(), [0.0, 7.388, 0.0].into(), 0.002564);
    let earth = MassiveParticle::new([1.0, 0.0, 0.0].into(), [0.0, 6.283, 0.0].into(), 0.003146);

    let asteroid = TestParticle::new([4.0, 0.0, 0.0].into(), [0.0, 7.0, 0.0].into());

    // Create particle lists
    let mut test_particles: Vec<TestParticle> = Vec::new();
    let mut massive_particles = Vec::new();

    // Push test data into particle lists
    massive_particles.push(sun);
    massive_particles.push(mercury);
    massive_particles.push(venus);
    massive_particles.push(earth);

    test_particles.push(asteroid);

    test_particles.shrink_to_fit();
    massive_particles.shrink_to_fit();

    // Initialise data read/write
    let mut wtr = WriterBuilder::new().from_path("data.csv")?;
    wtr.write_field("Time")?;
    wtr.write_field("")?;
    for _particle in test_particles.iter() {
        wtr.write_field("x")?;
        wtr.write_field("y")?;
        wtr.write_field("z")?;
        wtr.write_field("vx")?;
        wtr.write_field("vy")?;
        wtr.write_field("vz")?;
        wtr.write_field("")?;
    }
    for _particle in massive_particles.iter() {
        wtr.write_field("m")?;
        wtr.write_field("x")?;
        wtr.write_field("y")?;
        wtr.write_field("z")?;
        wtr.write_field("vx")?;
        wtr.write_field("vy")?;
        wtr.write_field("vz")?;
        wtr.write_field("")?;
    }
    wtr.write_record(None::<&[u8]>)?;

    // Compute System
    let mut current_time = 0.0;
    while current_time <= total_time {
        // Save data
        wtr.write_field(current_time.to_string())?;
        wtr.write_field("")?;
        for particle in test_particles.iter() {
            wtr.write_field(particle.position().x.to_string())?;
            wtr.write_field(particle.position().y.to_string())?;
            wtr.write_field(particle.position().z.to_string())?;
            wtr.write_field(particle.velocity().x.to_string())?;
            wtr.write_field(particle.velocity().y.to_string())?;
            wtr.write_field(particle.velocity().z.to_string())?;
            wtr.write_field("")?;
        }
        for particle in massive_particles.iter() {
            wtr.write_field(particle.mass().to_string())?;
            wtr.write_field(particle.position().x.to_string())?;
            wtr.write_field(particle.position().y.to_string())?;
            wtr.write_field(particle.position().z.to_string())?;
            wtr.write_field(particle.velocity().x.to_string())?;
            wtr.write_field(particle.velocity().y.to_string())?;
            wtr.write_field(particle.velocity().z.to_string())?;
            wtr.write_field("")?;
        }
        wtr.write_record(None::<&[u8]>)?;

        // Compute next load
        let next_output: f64 = current_time + output_interval;
        while current_time <= next_output {
            let mut test_accel = Vec::new();
            let mut mass_accel = Vec::new();

            for test_particle in test_particles.iter() {
                let mut acceleration = Vector3::<f64>::new(0.0, 0.0, 0.0);
                for massive_particle in massive_particles.iter() {
                    acceleration += test_particle.accel_towards(massive_particle);
                }
                test_accel.push(acceleration);
            }

            for massive_particle_a in massive_particles.iter() {
                let mut acceleration = Vector3::<f64>::new(0.0, 0.0, 0.0);
                for massive_particle_b in massive_particles.iter() {
                    if massive_particle_a != massive_particle_b {
                        acceleration += massive_particle_a.accel_towards(massive_particle_b);
                    }
                }
                mass_accel.push(acceleration);
            }

            for particle in test_particles.iter_mut().zip(test_accel.iter_mut()) {
                particle.0.update(particle.1, INTEGRATION_INTERVAL);
            }

            for particle in massive_particles.iter_mut().zip(mass_accel.iter_mut()) {
                particle.0.update(particle.1, INTEGRATION_INTERVAL);
            }

            current_time += INTEGRATION_INTERVAL;
        }
    }

    let elapsed = start.elapsed();
    println!("The program runtime was: {:?}", elapsed);
    wtr.flush()?;
    Ok(())
}
