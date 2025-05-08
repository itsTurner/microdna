import pysam
import argparse
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference',
                        help='Path to BAM file',
                        type=str,
                        required=True)
    parser.add_argument('-t', '--threshold',
                        help='Score threshold to include circle',
                        type=float,
                        default=1.0)
    parser.add_argument('-o', '--out_file',
                        help='Path to log file',
                        type=str)
    return parser.parse_args()


def is_softclipped(bam: pysam.AlignmentFile, read: pysam.AlignedSegment):
    return any(cigar_op[0] == pysam.CSOFT_CLIP for cigar_op in read.cigartuples)

def is_softclipped_with_match_ending(bam: pysam.AlignmentFile, read: pysam.AlignedSegment):
    return (read.cigartuples[-1][0] == pysam.CMATCH) and is_softclipped(bam, read)

def is_softclipped_with_match_beginning(bam: pysam.AlignmentFile, read: pysam.AlignedSegment):
    return (read.cigartuples[0][0] == pysam.CMATCH) and is_softclipped(bam, read)

def examine_softclipped_read(bam: pysam.AlignmentFile, read: pysam.AlignedSegment, score_threshold=1.0):
    blocks = read.get_blocks()
    end = blocks[0][0] + read.infer_query_length()
    potential_circles = []
    for block in blocks:
        if block[1] != end:
            for range_read in bam.fetch(reference=read.reference_name, start=block[1], end=end):
                if is_softclipped_with_match_beginning(bam, range_read):
                    overlap, depth, score = circle_metrics(start_read=read, end_read=range_read)
                    if score >= score_threshold:
                        potential_circles.append(range_read)
    return potential_circles

def report_potential_circles_for_completion(bam: pysam.AlignmentFile, read: pysam.AlignedSegment, potential_circles):
    report_string = f'{read.query_name}: {read.cigarstring} / {read.get_forward_sequence()} @ {read.get_blocks()[0][0]}\n'
    for range_read in potential_circles:
        start_str = read.get_forward_sequence()[:read.cigartuples[0][1]]
        end_str = range_read.get_forward_sequence()[read.cigartuples[-1][1]:]
        overlap, depth, score = circle_metrics(start_read=read, end_read=range_read)
        report_string += f'\tscore = {score:.2f} :: {range_read.cigarstring} / {range_read.get_forward_sequence()} / {range_read.get_blocks()[0][0]} / {end_str} <overlap> {start_str} = {overlap} / depth = {depth}\n'
    return report_string

def report_potential_circles_for_start_read(bam: pysam.AlignmentFile, read: pysam.AlignedSegment, potential_circles):
    print(f'{read.query_name}: {read.cigarstring} / {read.get_forward_sequence()} @ {read.get_blocks()}')
    for range_read in potential_circles:
        start_str = read.get_forward_sequence()[:read.cigartuples[0][1]]
        end_str = range_read.get_forward_sequence()[read.cigartuples[-1][1]:]
        overlap, depth, score = circle_metrics(start_read=read, end_read=range_read)
        print(f'\t{score:.2f} :: {range_read.cigarstring} / {range_read.get_forward_sequence()} / {range_read.get_blocks()} / {end_str} <> {start_str} = {overlap} / depth = {depth}')

def circle_metrics(start_read, end_read):
    overlap = amount_overlap(start_read=start_read, end_read=end_read)
    depth = depth_of_circle(start_read=start_read, end_read=end_read, overlap=overlap)
    score = overlap*depth/50
    return overlap, depth, score

def amount_overlap(start_read, end_read):
    start = start_read.get_forward_sequence()[:start_read.cigartuples[0][1]]
    end = end_read.get_forward_sequence()[end_read.cigartuples[-1][1]:]
    min_len = min(len(start), len(end))
    max_overlap = 0
    for i in range(min_len):
        if start[:i] == end[len(end)-i:]:
            max_overlap = i
    return max_overlap

def depth_of_circle(start_read, end_read, overlap):
    return start_read.infer_query_length() + end_read.cigartuples[-1][1] - overlap

def get_potential_circles(bam, reads, score_threshold=1.0):
    potential_circles = {}
    for read in reads:
        if is_softclipped_with_match_ending(bam, read):
            potential_circles_candidate = examine_softclipped_read(bam, read, score_threshold=1.0)
            if len(potential_circles_candidate) > 0:
                potential_circles[read] = potential_circles_candidate
    return potential_circles

def main():
    args = get_args()

    bam_file = args.reference

    if args.out_file and not os.path.exists(os.path.dirname(args.out_file)):
        raise ValueError('Out path given does not exist!')

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        final_report_string = ""
        potential_circles = get_potential_circles(bam, bam)
        for read, potential_circles_for_read in potential_circles.items():
            final_report_string += report_potential_circles_for_completion(bam, read, potential_circles_for_read)
        if args.out_file:
            with open(args.out_file, 'w') as f:
                f.write(final_report_string)
        else:
            print(final_report_string)

if __name__ == '__main__':
    main()