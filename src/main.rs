use std::path::PathBuf;

use clap::Parser;

use avio::{PixelFormat, VideoDecoder};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Video file to be used as input
    filename: PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let mut decoder = VideoDecoder::open(args.filename)
        .output_format(PixelFormat::Yuv420p)
        .build()?;

    while let Ok(Some(frame)) = decoder.decode_one() {}

    Ok(())
}
