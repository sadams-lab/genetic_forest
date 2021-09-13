pub fn mean(data: &Vec<&f64>) -> f64 {
    let sum: f64 = data.iter().filter_map(|value| 
        match **value {
            x if x >= 0. => Some(x),
            _ => None 
        }).sum();
    let count = data.len();
    match count {
        0 => return 0.,
        _ => {sum / count as f64}
    }

}

pub fn std_deviation(data: &Vec<&f64>) -> f64 {
    let data_mean = mean(data);
    let count = data.iter().filter_map(|value| 
        match **value {
            x if x >= 0. => Some(x),
            _ => None 
        }).count();
    match count {
        0 => return 0.,
        _ => {
            let variance = data.iter().filter_map(|value| 
                match **value {
                    x if x >= 0. => {
                        let diff = data_mean - x;
                        Some(diff * diff)
                    },
                    _ => None
            }).sum::<f64>() / count as f64;
            variance.sqrt()
        }
    }
}


