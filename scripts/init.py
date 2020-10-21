import os
import sys
import csv
import argparse
import collections
import ruamel.yaml
from pathlib import Path
from snakemake.utils import listfiles

def main():
    """Create a blank config file with the given name."""

    conda_prefix = get_conda_prefix()
    args = parse_args()
    project_fp = setup_project_folder(args)
    samplelists = write_samples_from_input(args, project_fp)
    write_config(args, conda_prefix, project_fp, samplelists)

def get_conda_prefix():
    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "Could not determine Conda prefix. Activate your ion-meta "
            "environment and try this command again.")
    return conda_prefix

def check_existing(path, force=False):
    if path.is_dir():
        raise SystemExit(
            "Error: specified file '{}' exists and is a directory".format(path))
    if path.is_file() and not force:
        raise SystemExit(
            "Error: specified file '{}' exists. Use --force to "
            "overwrite.".format(path))
    return path


def _write_samples_csv(samples, out):
    fieldnames = ["sample", "1", "2"]
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    for sample in samples.keys():
        writer.writerow({'sample':sample, **samples[sample]})

def build_sample_list(data_fp, format_str, output_file, is_single_end):

    data_fp = data_fp.resolve()
    fnames = [f.name for f in data_fp.iterdir() if f.is_file()]

    if not format_str:
        format_str = "{sample}.bam"

    samples = find_samples(data_fp, format_str)

    # Check for mate pairs if single end
    if not is_single_end:
        no_match = []
        for sample, reads in samples.items():
            if '2' not in reads:
                no_match.append(sample)
        if len(no_match) > 0:
            raise MissingMatePairError(
                "missing mate pairs for samples: {} ".format(
                    ", ".join(no_match)))

    if len(samples) == 0:
        raise SampleFormatError("no samples matching the given format found.")

    sys.stderr.write("Found {} samples in {}.\n".format(len(samples), data_fp))
    _write_samples_csv(samples, output_file)

def find_samples(data_fp, filename_fmt):
    """
    Build a list of samples from a data filepath and filename format.
    :param data_fp: a Path to data files
    :param filename_fmt: a string giving wildcards for {sample} and (optionally)
        {rp}, for read pair (e.g. R1 or R2).
    :returns: A dictionary of samples, with sample names as keys:
       Samples = {
         'sample1': {
           '1': 'path/to/sample1_R1.fastq.gz',
           '2': 'path/to/sample1_R2.fastq.gz', #optional
         }, ...
       }
    """
    files = list(listfiles(str(data_fp/filename_fmt)))
    Samples = {t[1]['sample']: {} for t in files}
    for f in files:
        fpath = f[0]
        wcards = f[1]
        rp = wcards.get('rp')
        if rp:
            if not rp in ['1', '2']:
                raise ValueError(
                    "'{rp}' should capture just '1' or '2' in filename, nothing else")
            Samples[wcards['sample']][rp] = fpath
        else:
            Samples[wcards['sample']]['1'] = fpath
    return Samples


def parse_args():
    description_str = (
        "Initialize a new ion-meta project in a given directory, creating "
        "a new config file and (optionally) a sample list.  The sample list "
        "source can be a folder of input files.")
    
    parser = argparse.ArgumentParser(
        "init", description=description_str)
    parser.add_argument(
        "project_fp", type=Path,
        help="project directory (will be created if it does not exist)")
    parser.add_argument(
        "-f", "--force", help="overwrite files if they already exist",
        action="store_true")
    parser.add_argument(
        "--output", help=(
            "name of config file (%(default)s)"),
        default="config.yaml", metavar="FILE")

    configs = parser.add_argument_group("config file options")
    configs.add_argument(
        "--defaults", type=argparse.FileType("r"), metavar="FILE",
        help="set of default values to use to populate config file")
    configs.add_argument(
        "--template", default="config/config.yaml", metavar="FILE",
        help="custom config file template, in YAML format", 
        type=argparse.FileType("r"))

    samplelist = parser.add_argument_group("sample list options",
            ("Options to automatically generate a sample list. --data_fp (for "
             "reading filenames from a specified folder) "
             ))
    
    samplelist.add_argument(
        "--data_fp", type=Path, metavar="PATH",
        help="path to folder containing .bam files")
    samplelist.add_argument(
        "--format", type=str, metavar="STR",
        help="filename format for --data_fp (default: IonXpress_{sample}.bam)")
    samplelist.add_argument(
        "--single_end", action="store_true",
        help="bam files are in single-end, not paired-end")

    args = parser.parse_args()

    return args

def setup_project_folder(args):
    # Create project folder if it doesn't exist
    project_fp_exists = False
    project_fp = args.project_fp
    try:
        project_fp = args.project_fp.resolve()
        project_fp_exists = project_fp.exists()
    except FileNotFoundError:
        pass
    if not project_fp_exists:
        sys.stderr.write(
            "Creating project folder at {}...\n".format(args.project_fp))
        project_fp.mkdir(parents=True, exist_ok=True)
    return project_fp

def write_samples_from_input(args, project_fp):
    """Write sample list CSV from existing files."""
    samplelist_file = check_existing(project_fp/"samples.csv", args.force)
    if args.data_fp:
            try:
                with samplelist_file.open("w") as out:
                    build_sample_list(
                        data_fp = args.data_fp,
                        format_str = args.format,
                        output_file = out,
                        is_single_end = args.single_end)
            except SampleFormatError as e:
                raise SystemExit(
                    "Error: could not create sample list. Specify correct sample filename"
                    " format using --format.\n  Reason: {}".format(e))
            except MissingMatePairError as e:
                raise SystemExit(
                    "Error: assuming paired-end reads, but could not find mates. Specify "
                    "--single_end if not paired-end, or provide sample name format "
                    "using --format."
                    "\n  Reason: {}".format(e))
            sys.stderr.write(
                "New sample list written to {}\n".format(samplelist_file))

    samplelists = {["paired", "unpaired"][args.single_end]: samplelist_file}

    return samplelists

def config_new(
        conda_fp, project_fp,
        template=None):
    if template:
        config = template.read()
    return config.format(
        CONDA_FP=conda_fp,
        PROJECT_FP=project_fp)

def _update_dict(target, new):
    for k, v in new.items():
        if isinstance(v, collections.Mapping):
            # We could use .get() here but ruamel.yaml's weird Mapping
            # subclass outputs errors to stdout if the key doesn't exist
            if k in target:
                target[k] = _update_dict(target[k], v)
            else:
                target[k] = _update_dict({}, v)
        else:
            target[k] = v
    return target

def _update_dict_strict(target, new):
    for k, v in new.items():
        if isinstance(v, collections.Mapping) and k in target.keys():
            target[k] = _update_dict_strict(target.get(k, {}), v)
        elif k in target.keys():
            target[k] = v
        else:
            sys.stderr.write("Key '%s' not found in target, skipping\n" % k)
            continue
    return target

def config_update(config_str, new, strict=False):
    config = ruamel.yaml.round_trip_load(config_str)
    if strict:
        config = _update_dict_strict(config, new)
    else:
        config = _update_dict(config, new)
    return config

def config_dump(config, out=sys.stdout):
    if isinstance(config, collections.Mapping):
        ruamel.yaml.round_trip_dump(config, out)
    else:
        out.write(config)

def write_config(args, conda_prefix, project_fp, samplelists):
    multiple_configs = len(samplelists.values()) != 1
    for layout in samplelists.keys(): # one paired, one unpaired, or one of each
        if multiple_configs:
            config_file = check_existing(project_fp/Path(layout+"_"+str(Path(project_fp/args.output).name)), args.force)
        else:
            config_file = check_existing(project_fp/args.output, args.force)
        cfg = config_new(
            conda_fp=conda_prefix,
            project_fp=project_fp,
            template=args.template)
        defaults = {}
        if args.defaults:
            defaults = ruamel.yaml.safe_load(args.defaults)
        # Override loaded config defaults (if any) for a few specific items.
        paired = layout == "paired"
        defaults["paired_end"] = paired
        #print(type(samplelists)) # <class 'dict'>
        #print(samplelists[layout])
        #print(samplelists[layout].name)
        #defaults["samplelist_fp"] = samplelists[layout].name

        cfg = config_update(cfg, defaults)
        with config_file.open('w') as out:
            config_dump(cfg, out)
        sys.stderr.write("New config file written to {}\n".format(config_file))
    if multiple_configs:
        raise SystemExit("Found both paired and unpaired reads. Wrote two sample lists "
                        "and config files, with '_paired' or '_single' appended.")

class MissingMatePairError(Exception):
    pass

class SampleFormatError(Exception):
    pass

if __name__ == '__main__':
    main()
