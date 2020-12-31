pub struct Variant {
    id: String,
    max_importance: f32
}

impl Variant {
    pub fn new(id: String) -> Self {
        return Variant {
            id: id, 
            max_importance = 0.
        }
    }
    pub fn update_importance(&self, imp: f32) {
        &self.max_importance = imp
    }
}