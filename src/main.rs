use std::error::Error;
use std::time::Instant;
use N_orbital::*;
// use nalgebra::Vector3;
use csv::WriterBuilder;
use itertools::Itertools;

const SAVE_RESULTS: bool = true; // This is the wrong way to do this, but that's okay

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
    if SAVE_RESULTS {
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
    }

    // Compute System
    let mut current_time = 0.0;
    while current_time <= total_time {
        // Save data
        if SAVE_RESULTS {
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
        }

        // Compute next load
        let next_output: f64 = current_time + output_interval;
        while current_time <= next_output {

            for (test, massive) in test_particles
                .iter()
                .cartesian_product(massive_particles.iter())
            {
                test_calc(test, massive)?;
            }

            for pair in massive_particles.iter().combinations(2) {
                massive_calc(pair[0], pair[1])?;
            }

            for particle in test_particles.iter_mut() {
                particle.update(INTEGRATION_INTERVAL)?;
            }

            for particle in massive_particles.iter_mut() {
                particle.update(INTEGRATION_INTERVAL)?;
            }

            current_time += INTEGRATION_INTERVAL;
        }
    }

    let elapsed = start.elapsed();
    println!("The program runtime was: {:?}", elapsed);
    if SAVE_RESULTS {
        wtr.flush()?;
    }
    Ok(())
}
