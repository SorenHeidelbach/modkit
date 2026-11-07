use anyhow::Result as AnyhowResult;
use clap::Subcommand;
use mod_kit::pileup::subcommand::ModBamPileup;

#[derive(Subcommand)]
pub enum Commands {
    /// Tabulates base modification calls across genomic positions. This
    /// command produces a bedMethyl formatted file. Schema and description
    /// of fields can be found in the README.
    Pileup(ModBamPileup),
}

impl Commands {
    pub fn run(&self) -> AnyhowResult<()> {
        match self {
            Self::Pileup(x) => x.run(),
        }
    }
}
