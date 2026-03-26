#!/usr/bin/env python3
"""Generate pairwise combinatorial test configs for VarDict CLI options.

Uses allpairspy to create a pairwise covering array from the PICT model.
Install: pip install allpairspy

Usage:
    python tests/generate_pairwise_configs.py > tests/pairwise_configs.tsv
    python tests/generate_pairwise_configs.py --shell > tests/pairwise_option_parity.sh

The output TSV can be fed into the option parity test runner.
"""
import sys
from itertools import product

try:
    from allpairspy import AllPairs
except ImportError:
    print("ERROR: Install allpairspy first: pip install allpairspy", file=sys.stderr)
    sys.exit(1)

# Parameter definitions: name -> (flag, values_dict)
# values_dict maps symbolic values to CLI flag strings
PARAMETERS = {
    "pileup":         ("-p",  {"off": "",    "on": "-p"}),
    "freq":           ("-f",  {"0.01": "-f 0.01", "0.001": "-f 0.001", "0.05": "-f 0.05"}),
    "minr":           ("-r",  {"2": "-r 2",  "1": "-r 1",  "5": "-r 5"}),
    "qual":           ("-q",  {"22.5": "",   "0": "-q 0",  "30": "-q 30"}),
    "mapq":           ("-Q",  {"0": "",      "10": "-Q 10", "30": "-Q 30"}),
    "realign":        ("-k",  {"1": "",      "0": "-k 0"}),
    "nosv":           ("-U",  {"off": "",    "on": "-U"}),
    "fisher":         ("--fisher", {"off": "", "on": "--fisher"}),
    "vext":           ("-X",  {"2": "",      "0": "-X 0",  "5": "-X 5"}),
    "mismatch":       ("-m",  {"8": "",      "3": "-m 3",  "15": "-m 15"}),
    "filter":         ("-F",  {"0x504": "",  "0x700": "-F 0x700", "0x100": "-F 0x100"}),
    "move_3prime":    ("-3",  {"off": "",    "on": "-3"}),
    "read_pos":       ("-P",  {"5.0": "",    "0": "-P 0",  "10": "-P 10"}),
    "qratio":         ("-o",  {"1.5": "",    "0.5": "-o 0.5", "2.0": "-o 2.0"}),
    "chimeric":       ("--chimeric", {"off": "", "on": "--chimeric"}),
    "deldupvar":      ("--deldupvar", {"off": "", "on": "--deldupvar"}),
    "include_n":      ("-K",  {"off": "",    "on": "-K"}),
    "extend":         ("-x",  {"0": "",      "150": "-x 150"}),
    "ref_extension":  ("-Y",  {"1200": "",   "600": "-Y 600", "2000": "-Y 2000"}),
    "insert_size":    ("-w",  {"300": "",    "200": "-w 200"}),
    "insert_std":     ("-W",  {"100": "",    "50": "-W 50"}),
    "trim_bases":     ("-T",  {"0": "",      "130": "-T 130"}),
    "min_match":      ("-M",  {"0": "",      "25": "-M 25"}),
    "unique_mode":    ("-u",  {"off": "",    "u": "-u", "UN": "-UN"}),
    "insert_std_amt": ("-A",  {"4": "",      "2": "-A 2", "6": "-A 6"}),
    "sv_min_len":     ("-L",  {"1000": "",   "500": "-L 500", "2000": "-L 2000"}),
    "min_bias_reads": ("-B",  {"2": "",      "1": "-B 1", "4": "-B 4"}),
}

PARAM_NAMES = list(PARAMETERS.keys())
PARAM_VALUES = [list(PARAMETERS[name][1].keys()) for name in PARAM_NAMES]


def is_valid_combination(row):
    """Filter function for constraints matching PICT model."""
    if len(row) < 2:
        return True
    
    vals = {}
    for i, v in enumerate(row):
        vals[PARAM_NAMES[i]] = v
    
    # pileup=on requires freq <= 0.01
    if vals.get("pileup") == "on" and vals.get("freq") == "0.05":
        return False
    
    # nosv=on makes SV params irrelevant (force defaults)
    if vals.get("nosv") == "on":
        if vals.get("insert_size") not in (None, "300"):
            return False
        if vals.get("insert_std") not in (None, "100"):
            return False
        if vals.get("sv_min_len") not in (None, "1000"):
            return False
        if vals.get("insert_std_amt") not in (None, "4"):
            return False
    
    return True


def generate_configs():
    """Generate pairwise covering array."""
    pairs = AllPairs(PARAM_VALUES, filter_func=is_valid_combination)
    configs = []
    for i, row in enumerate(pairs):
        config = {}
        for j, val in enumerate(row):
            config[PARAM_NAMES[j]] = val
        configs.append(config)
    return configs


def config_to_cli_flags(config):
    """Convert a config dict to CLI flag string."""
    flags = []
    for name, value in config.items():
        _, values_dict = PARAMETERS[name]
        flag_str = values_dict[value]
        if flag_str:
            flags.append(flag_str)
    return " ".join(flags) if flags else "(default)"


def config_to_label(config, index):
    """Generate a human-readable label for a config."""
    non_default = []
    for name, value in config.items():
        values = list(PARAMETERS[name][1].keys())
        if value != values[0]:  # Not the default value
            non_default.append(f"{name}={value}")
    if not non_default:
        return f"PW-{index:03d}-default"
    return f"PW-{index:03d}-{'_'.join(non_default[:3])}"


def main():
    shell_mode = "--shell" in sys.argv
    configs = generate_configs()
    
    if shell_mode:
        print("#!/usr/bin/env bash")
        print("# Auto-generated pairwise configs from generate_pairwise_configs.py")
        print(f"# {len(configs)} configs covering all 2-way parameter interactions")
        print(f"# Parameters: {len(PARAM_NAMES)}")
        print()
        print("PAIRWISE_CONFIGS=(")
        for i, config in enumerate(configs):
            label = config_to_label(config, i)
            flags = config_to_cli_flags(config)
            print(f'  "{label}|{flags}"')
        print(")")
    else:
        # TSV header
        print("\t".join(["index", "label", "cli_flags"] + PARAM_NAMES))
        for i, config in enumerate(configs):
            label = config_to_label(config, i)
            flags = config_to_cli_flags(config)
            values = [config[name] for name in PARAM_NAMES]
            print("\t".join([str(i), label, flags] + values))
    
    print(f"\n# Total configs: {len(configs)}", file=sys.stderr)
    print(f"# Parameters: {len(PARAM_NAMES)}", file=sys.stderr)
    print(f"# This covers all 2-way interactions between {len(PARAM_NAMES)} parameters", file=sys.stderr)


if __name__ == "__main__":
    main()
