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


