# Documentation Index

This directory contains detailed documentation for the L-BFGS-B Rust port.

## Files

### [EXAMPLES.md](EXAMPLES.md)
Comprehensive guide to the three example drivers:
- `driver1.rs` - Basic usage
- `driver2.rs` - Custom stopping criteria
- `driver3.rs` - Time-controlled optimization

Includes API examples, callback structure documentation, and usage instructions.

### [DRIVER_PORT_SUMMARY.md](DRIVER_PORT_SUMMARY.md)
Technical summary of the driver2 and driver3 port from C to Rust:
- Implementation details
- API design rationale
- Testing and verification results
- Comparison with C reference implementation

### [COMPARISON.md](COMPARISON.md)
Detailed comparison between the C reference implementation and the Rust port:
- File-by-file mapping
- Function equivalence tables
- Implementation differences
- Porting notes

### [PORTING_STATUS.md](PORTING_STATUS.md)
Status document tracking the porting progress:
- Completed components
- Known issues and fixes
- Testing status
- Implementation notes

### [SUMMARY.md](SUMMARY.md)
Executive summary of the entire porting project:
- Overview of changes
- Key accomplishments
- Verification results
- Next steps and recommendations

## Quick Start

For usage examples, start with [EXAMPLES.md](EXAMPLES.md).

For understanding the port's fidelity to the original, see [COMPARISON.md](COMPARISON.md).

For technical details about the driver ports, see [DRIVER_PORT_SUMMARY.md](DRIVER_PORT_SUMMARY.md).