#! /usr/bin/env python
__author__ = 'BinglanLi'

import allel
import gzip
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.parse
import urllib.request

import vcf_preprocess_exceptions as Exceptions


# chromosome names
_chr_invalid = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
                "18", "19", "20", "21", "22", "X", "Y", "M", "MT", "chrMT"]
_chr_valid = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
              "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
              "chr22", "chrX", "chrY", "chrM", "chrM", "chrM"]
_chr_valid_sorter = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                     "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                     "chr22", "chrX", "chrY", "chrM"]
# check if two chr arrays are of the same length
if len(_chr_invalid) != len(_chr_valid):
    print("Error in internal chromosome mapping arrays")
    sys.exit(1)


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def bgzip_file(bgzip_path, vcf_path):
    """
    bgzip the file
    """

    try:
        print("Bgzipping", vcf_path)
        subprocess.run([bgzip_path, '-f', vcf_path], check=True)
    except Exception as e:
        print('Failed to bgzip %s' % vcf_path)
        # comment out this traceback function as subprocess(check = true) should report detailed errors
        # traceback.print_exception(type(e), e, e.__traceback__)
        sys.exit(1)


def bgzipped_vcf(bgzip_path, file):
    """
    make sure file is bgzipped
    """
    if not is_gz_file(file):
        print("VCF file is not bgzipped!")
        bgzip_file(bgzip_path, file)
        file = file + '.gz'
        if os.path.exists(file + '.tbi'):
            print("Removing pre-existing .tbi")
            os.remove(file + '.tbi')
    return file


def byte_decoder(a):
    """ Decode byte data into utf-8 """
    return a.decode("utf-8")


def download_from_url(url, download_to_dir, save_to_file=None, force_update=False):
    """download from an url"""

    remote_basename = os.path.basename(urllib.parse.urlparse(url).path)
    if remote_basename:
        local_path = os.path.join(download_to_dir, remote_basename) if not save_to_file else save_to_file
        if os.path.exists(local_path):
            if force_update:
                os.remove(local_path)
            else:
                return local_path

        # Download to a temp file. If a download succeeds, rename the temp file.
        # If a download fails, the function will throw an exception. The temp file will be removed.
        with tempfile.TemporaryDirectory(dir=download_to_dir) as temp_dir:
            temp_download_path = os.path.join(temp_dir, 'temp_' + remote_basename)
            with urllib.request.urlopen(url) as response:
                with open(temp_download_path, 'wb') as out_file:
                    print('Downloading from \"%s\"\n\tto \"%s\"' % (url, local_path))
                    shutil.copyfileobj(response, out_file)
            os.rename(temp_download_path, local_path)
        return local_path
    else:
        raise Exceptions.InvalidURL(url)


def get_default_grch38_ref_fasta_and_index(download_to_dir, force_update=False):
    """download the human reference genome sequence GRCh38/hg38"""

    ref_file = os.path.join(download_to_dir, 'reference.fasta.bgz')
    if os.path.exists(ref_file) and not force_update:
        return ref_file

    tar_file = download_from_url('https://zenodo.org/record/5572839/files/GRCh38_reference_fasta.tar?download=1',
                                 download_to_dir, None, force_update)
    with tarfile.open(tar_file, 'r') as tar:
        tar.extractall(path=download_to_dir)
    os.remove(tar_file)

    return ref_file


def tabix_index_vcf(tabix_path, vcf_path):
    """
    index the input vcf using tabix, and the output index file will be written to the working directory

    tabix commands are exclusively "tabix -p vcf <input_vcf>", which generates an index file (.tbi)
        for an input file (<input_file>) whose file type is specified by "-p vcf".
    .tbi will be output to the current working directory by default.
    """

    try:
        print('Generating index (' + vcf_path + '.tbi)')
        subprocess.run([tabix_path, '-p', 'vcf', vcf_path], check=True)
    except Exception as e:
        print('Failed to index %s' % vcf_path)
        # comment out this traceback function as subprocess(check = true) should report detailed errors
        # traceback.print_exception(type(e), e, e.__traceback__)
        sys.exit(1)


def run_bcftools(list_bcftools_command, show_msg=None):
    """
    run the bcftools following the commands stored in the list_bcftools_command

    "bcftools <common_options> <input_vcf>".
    "-Oz" (capitalized letter O) specifies the output type as compressed VCF (z). "-o" writes to a file rather than to
    default standard output.
    "--no-version" will cease appending version and command line information to the output VCF header.
    "-s sample_ID(s)" comma-separated list of samples to include or exclude if prefixed with "^".
    "-r chr|chr:pos|chr:beg-end|chr:beg-[,…]" extracts comma-separated list of regions
    """

    print("%s" % show_msg) if show_msg else print("Running [ %s ]" % (' '.join(list_bcftools_command)))

    p = subprocess.run(list_bcftools_command, stderr=subprocess.PIPE)
    if p.returncode != 0:
        print(p.stderr.decode("utf-8"))
        sys.exit(1)

def obtain_vcf_sample_list(bcftools_path, path_to_vcf):
    """
    obtain a list of samples from the input VCF

    "bcftools query <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-l" list sample names and exit.
    Samples are delimited by '\\\\n' and the last line ends as 'last_sample_ID\\\\n\\\\n'.
    """

    output = subprocess.check_output([bcftools_path, 'query', '-l', path_to_vcf], universal_newlines=True)
    vcf_sample_list = output.split('\n')[:-1]  # remove the black line at the end
    return vcf_sample_list


def remove_vcf_and_index(path_to_vcf):
    """remove the compressed vcf as well as the index file"""

    try:
        os.remove(path_to_vcf)
        os.remove(path_to_vcf + '.tbi')
        print("Removed intermediate files:\n\t%s\n\t%s" % (path_to_vcf, path_to_vcf + '.tbi'))
    except OSError as error_remove_tmp:
        print("Error: %s : %s" % (path_to_vcf, error_remove_tmp.strerror))


def _get_vcf_pos_min_max(positions, flanking_bp=100):
    """ given input positions, return "<min_pos>-<max_pos>"  """
    return '-'.join([str(min(positions) - flanking_bp), str(max(positions) + flanking_bp)])


def _is_valid_chr(input_vcf):
    chr_status = 'none'
    with gzip.open(input_vcf) as f:
        for line in f:
            try:
                line = byte_decoder(line)
            except:
                line = line
            if line[0] != '#':
                fields = line.rstrip().split()
                if fields[0] in _chr_valid:
                    chr_status = 'true'
                    break
                elif fields[0] in _chr_invalid:
                    chr_status = 'false'
                    break
                else:
                    break
    return chr_status


def extract_regions_from_single_file(bcftools_path, tabix_path, input_vcf, pgx_vcf, output_dir, output_prefix,
                                     sample_list):
    """
    Rename chromosomes in input vcf according to a chr-renaming mapping file.
    Extract pgx regions from input_vcf based on the pgx_vcf.

    "bcftools annotate <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "--rename-chrs" renames chromosomes according to the map in file_rename_chrs.
    """

    # output path
    path_output = os.path.join(output_dir, 'PharmCAT_preprocess_' + output_prefix + '.pgx_regions.vcf.gz')

    print("")
    print("Processing", input_vcf)
    # create index if not already existed
    if not os.path.exists(input_vcf + '.tbi'):
        tabix_index_vcf(tabix_path, input_vcf)

    # obtain PGx regions to be extracted
    df_ref_pgx_pos = allel.vcf_to_dataframe(pgx_vcf)
    df_ref_pgx_pos['CHROM'] = df_ref_pgx_pos['CHROM'].astype("category")
    df_ref_pgx_pos['CHROM'].cat.set_categories(_chr_valid_sorter, inplace=True)
    ref_pgx_regions = df_ref_pgx_pos.groupby(['CHROM'])['POS'].agg(_get_vcf_pos_min_max).reset_index()
    ref_pgx_regions.dropna(axis=0, subset=['POS'], how='any', inplace=True)
    # add a special case for 'chrMT'
    idx_chrM = ref_pgx_regions.index[ref_pgx_regions['CHROM'] == 'chrM']
    ref_pgx_regions = ref_pgx_regions.append(
        ref_pgx_regions.loc[idx_chrM].assign(**{'CHROM': 'chrMT'}), ignore_index=True)

    # generate a temp dir to extract pgx regions and, if necessary, rename chr
    with tempfile.TemporaryDirectory(suffix='extract_pgx_regions', dir=output_dir) as temp_dir:
        # generate temp file of sample list
        file_sample_list = os.path.join(temp_dir, 'sample_list.txt')
        with open(file_sample_list, 'w+') as f:
            for single_sample in sample_list:
                f.write("%s\n" % single_sample)

        # create a temporary chromosome mapping file
        file_chr_rename = os.path.join(temp_dir, 'rename_chr.txt')
        with open(file_chr_rename, 'w+') as f:
            for i in range(len(_chr_invalid)):
                f.write(_chr_invalid[i] + "\t" + _chr_valid[i] + "\n")

        # validate chromosome formats
        if _is_valid_chr(input_vcf) == 'true':
            # format the pgx regions to be extracted
            ref_pgx_regions = ",".join(ref_pgx_regions.apply(lambda row: ':'.join(row.values.astype(str)), axis=1))
        elif _is_valid_chr(input_vcf) == 'false':
            # format pgx regions
            ref_pgx_regions = ",".join(
                ref_pgx_regions.apply(lambda row: ':'.join(row.values.astype(str)), axis=1).replace({'chr': ''},
                                                                                                    regex=True))
        else:
            print("The CHROM column does not comply with either 'chr##' or '##' format.")
            sys.exit(1)

        # extract pgx regions and modify chromosome names if necessary
        bcftools_command = [bcftools_path, 'annotate', '--no-version', '-S', file_sample_list,
                            '--rename-chrs', file_chr_rename, '-r', ref_pgx_regions,
                            '-Oz', '-o', path_output, input_vcf]
        run_bcftools(bcftools_command,
                     show_msg='Extracting PGx regions and modifying chromosome names for %s.' % input_vcf)

    # index the output PGx file
    tabix_index_vcf(tabix_path, path_output)

    return path_output


def extract_regions_from_multiple_files(bcftools_path, tabix_path, bgzip_path, input_list, pgx_vcf,
        output_dir, output_prefix, sample_list):
    """
    iterate through the list of input files
    """

    path_output = os.path.join(output_dir, 'PharmCAT_preprocess_' + output_prefix + '.pgx_regions.vcf.gz')

    with tempfile.TemporaryDirectory(suffix='concat_input_vcfs', dir=output_dir) as temp_dir:
        # process the vcfs in the input list one by one
        preprocessed_file_list = []
        i = 1
        with open(input_list, 'r') as file:
            for line in file:
                line = line.strip()
                if os.path.isfile(line):
                    print("")
                    print("Processing", line)
                    if not os.path.exists(line):
                        print("Cannot find", line)
                        continue
                    line = bgzipped_vcf(bgzip_path, line)

                    temp_output_prefix = output_prefix + '_' + i
                    single_file = extract_regions_from_single_file(bcftools_path, tabix_path, line, pgx_vcf,
                                                                   temp_dir, temp_output_prefix, sample_list)
                    preprocessed_file_list.append(single_file)
                    i += 1
                else:
                    print("Warning: Skip %s because the file does not exist." % line)

        # generate a temporary list of files to be concatenated
        temp_file_list = os.path.join(temp_dir, "temp_file_list.txt")
        with open(temp_file_list, 'w+') as f:
            for j in range(len(preprocessed_file_list)):
                f.write(preprocessed_file_list[j] + "\n")

        # concatenate vcfs
        bcftools_command = [bcftools_path, 'concat', '--no-version', '-a', '-f', temp_file_list, '-Oz', '-o',
                            path_output]
        run_bcftools(bcftools_command, show_msg='Concatenating chromosome VCFs.')

    # index the concatenated VCF
    tabix_index_vcf(tabix_path, path_output)
    return path_output


def normalize_vcf(bcftools_path, tabix_path, input_vcf, ref_seq):
    """
    Normalize the input VCF against the human reference genome sequence GRCh38/hg38

    "bcftools norm <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-m +|-" joins biallelic sites into multiallelic records (+)
        and convert multiallelic records into uniallelic format (-).
    "-f <ref_seq_fasta>" reference sequence. Supplying this option turns on left-alignment and normalization.
    "-c ws" when incorrect or missing REF allele is encountered, warn (w) and set/fix(s) bad sites.  's' will swap
    alleles and update GT and AC acounts. Importantly, s will NOT fix strand issues in a VCF.
    """

    path_output = os.path.splitext(os.path.splitext(input_vcf)[0])[0] + '.normalized.vcf.gz'

    bcftools_command = [bcftools_path, 'norm', '--no-version', '-m-', '-c', 'ws', '-Oz', '-o',
                        path_output, '-f', ref_seq, input_vcf]
    run_bcftools(bcftools_command, show_msg='Normalizing VCF')
    tabix_index_vcf(tabix_path, path_output)

    return path_output


def filter_pgx_variants(bcftools_path, tabix_path, input_vcf, ref_seq, pgx_vcf, output_dir, output_prefix):
    """
    Extract specific pgx positions that are present in the reference PGx VCF
    Generate a report of PGx positions that are missing in the input VCF

    "bcftools isec <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-c none" only records with the same CHR, POS, REF and ALT are considered identical
    "-C" outputs positions present only in the first file but missing in the others.
    "-w" lists input files to output given as 1-based indices.
        "-w1" extracts and writes records only present in the first file (the reference PGx positions).
    """

    path_output = os.path.splitext(os.path.splitext(input_vcf)[0])[0] + '.multiallelic.vcf.gz'

    with tempfile.TemporaryDirectory(suffix='extract_pgx_variants', dir=output_dir) as temp_dir:
        # convert reference PGx variants to the uniallelic format
        ref_pgx_uniallelic = os.path.join(temp_dir, 'temp_reference_pgx_variants_sorted_uniallelic.vcf.gz')
        bcftools_command = [bcftools_path, 'norm', '--no-version', '-m-', '-c', 'ws', '-f', ref_seq,
                            '-Oz', '-o', ref_pgx_uniallelic, pgx_vcf]
        run_bcftools(bcftools_command, show_msg='Preparing the reference PGx VCF')
        tabix_index_vcf(tabix_path, ref_pgx_uniallelic)

        # extract only the required PGx positions
        input_pgx_variants_only = os.path.join(temp_dir, 'temp_input_pgx_variants_only.vcf.gz')
        # output the exact matching variants (both position and alleles) in the second input file
        bcftools_command = [bcftools_path, 'isec', '--no-version', '-c', 'none', '-n=2', '-w2',
                            '-Oz', '-o', input_pgx_variants_only,
                            ref_pgx_uniallelic, input_vcf]
        run_bcftools(bcftools_command, show_msg='Retaining only PGx positions')
        tabix_index_vcf(tabix_path, input_pgx_variants_only)

        # merging the filtered input with the reference PGx positions for the next step
        merge_vcf = os.path.join(temp_dir, 'temp_merged.vcf.gz')
        bcftools_command = [bcftools_path, 'merge', '--no-version', '-m', 'both',
                            '-Oz', '-o', merge_vcf, input_pgx_variants_only, ref_pgx_uniallelic]
        run_bcftools(bcftools_command,
                     show_msg='Trimming file')
        tabix_index_vcf(tabix_path, merge_vcf)

        # enforce the output to comply with PharmCAT format
        multiallelic_vcf = os.path.join(temp_dir, 'temp_multiallelic.vcf.gz')
        bcftools_command = [bcftools_path, 'norm', '--no-version', '-m+', '-N',
                            '-Oz', '-o', multiallelic_vcf, merge_vcf]
        run_bcftools(bcftools_command,
                     show_msg='Enforcing the variant representation per PharmCAT')
        tabix_index_vcf(tabix_path, multiallelic_vcf)

        # remove the artificial PharmCAT sample
        bcftools_command = [bcftools_path, 'view', '--no-version', '-s', '^PharmCAT',
                            '-Oz', '-o', path_output, multiallelic_vcf]
        run_bcftools(bcftools_command, show_msg='Trimming file')
        tabix_index_vcf(tabix_path, path_output)

        # report missing positions in the input VCF
        report_missing_variants = os.path.join(output_dir, output_prefix + '.missing_pgx_var.vcf.gz')
        bcftools_command = [bcftools_path, 'isec', '--no-version', '-c', 'none', '-w1',
                            '-Oz', '-o', report_missing_variants, '-C', ref_pgx_uniallelic, input_vcf]
        run_bcftools(bcftools_command,
                     show_msg='Generating a report of missing PGx allele defining positions')

    return path_output


def output_pharmcat_ready_vcf(bcftools_path, input_vcf, output_dir, output_prefix, sample_list):
    """
    iteratively write to a PharmCAT-ready VCF for each sample

    "bcftools view <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-U" exclude sites without a called genotype, i.e., GT = './.'
    """

    for single_sample in sample_list:
        output_file_name = os.path.join(output_dir, output_prefix + '.' + single_sample + '.vcf')
        bcftools_command = [bcftools_path, 'view', '--no-version', '-U', '-Ov',
                            '-o', output_file_name, '-s', single_sample, input_vcf]
        run_bcftools(bcftools_command,
                     show_msg='Generating a PharmCAT-ready VCF for ' + single_sample)
