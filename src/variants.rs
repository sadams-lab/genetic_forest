// Copyright 2020 Solomon M. Adams, PharmD, PhD
// Licensed under the MIT license

//! Manage indices and meta-data about variants
//! used in the random forest algorithm

pub struct Variant {
    pub id: String,
    pub max_importance: f64
}

impl Variant {
    pub fn new(id: String) -> Self {
        return Variant {
            id: id, 
            max_importance: 0.
        }
    }
    pub fn set_importance(&mut self, imp: f64) {
        self.max_importance = imp
    }
}