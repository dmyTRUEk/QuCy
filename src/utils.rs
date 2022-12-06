/// Utils

use std::fmt::Display;


pub fn exit_with_msg(msg: &str) -> ! {
    println!("{msg}");
    std::process::exit(0);
}

pub trait ExtensionUnwrapOrExit<T> {
    fn unwrap_or_exit_with_msg(self, msg: &str) -> T;
}
impl<T> ExtensionUnwrapOrExit<T> for Option<T> {
    fn unwrap_or_exit_with_msg(self, msg: &str) -> T {
        self.unwrap_or_else(|| exit_with_msg(&format!("{msg}.")))
    }
}
impl<T, E: Display> ExtensionUnwrapOrExit<T> for Result<T, E> {
    fn unwrap_or_exit_with_msg(self, msg: &str) -> T {
        self.unwrap_or_else(|e| exit_with_msg(&format!("{msg}: {e}.")))
    }
}

