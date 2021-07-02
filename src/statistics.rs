// Statistical measures
use crate::z_table;

pub fn mean(data: &Vec<&f64>) -> f64 {
    let sum: f64 = data.iter().map(|value| **value).sum();
    let count = data.len();
    match count {
        0 => return 0.,
        _ => {sum / count as f64}
    }

}

pub fn std_deviation(data: &Vec<&f64>) -> f64 {
    let data_mean = mean(data);
    let count = data.len();
    match count {
        0 => return 0.,
        _ => {
            let variance = data.iter().map(|value| {
                let diff = data_mean - (**value as f64);
                diff * diff
            }).sum::<f64>() / count as f64;
            variance.sqrt()
        }
    }
}

pub fn get_cutoff(sd: &f64, mean: &f64, p: &f64) -> f64 {
    let mut p_index: usize = 0;
    if *p > 0.5 {
        return 0.
    }
    while p_index < (z_table::MAX_Z * 100.0) as usize {
        if (z_table::TABLE[p_index] as f64 * 100.0).round() == ((1.0 - p) * 100.0).round() {
            return (p_index as f64 / 100.0) * sd + mean
        } 
        p_index += 1
    }
    4.0 * sd + mean
    
}

