use anyhow::{anyhow, Result as AnyhowResult};
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::process::Output;

pub fn run_modkit(args: &[&str]) -> AnyhowResult<Output> {
    let exe = Path::new(env!("CARGO_BIN_EXE_modkit"));
    assert!(exe.exists());

    // eprintln!("running modkit with args {:?}", args.join(" "));
    let output = std::process::Command::new(exe)
        .args(args)
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()?
        .wait_with_output()?;
    if output.status.success() {
        Ok(output)
    } else {
        Err(anyhow!("failed to run {:?}", args.join(" ")))
    }
}

pub fn check_against_expected_text_file(output_fp: &str, expected_fp: &str) {
    assert_ne!(output_fp, expected_fp, "cannot check a file against itself");
    let test = {
        let mut fh = File::open(output_fp).unwrap();
        let mut buff = String::new();
        fh.read_to_string(&mut buff).unwrap();
        buff
    };
    let expected = {
        // this file was hand-checked for correctness or should be equivalent
        // due to orthogonal process
        let mut fh = File::open(expected_fp).unwrap();
        let mut buff = String::new();
        fh.read_to_string(&mut buff).unwrap();
        buff
    };

    similar_asserts::assert_eq!(
        test,
        expected,
        "{output_fp} is not the same as {expected_fp}"
    );
}
