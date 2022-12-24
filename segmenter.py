import argparse
import json
from Bio import pairwise2
import editdistance

GAP_CHAR = "\u2581"


def get_text(jsons, key, keep_lines=False):
    words = []
    for line in jsons:
        lw = line[key].split()
        if keep_lines:
            words.append(lw)
        else:
            words += lw

    return words


def match_fn(a, b):
    ml = max(len(a), len(b))
    return 1 - editdistance.eval(a, b) / ml


def match_fn_mwer(a, b):
    return a == b


def segment(
    unsegmented,
    segmented,
    json_unsegmented_transcript_key,
    json_segmented_transcript_key,
    mwer_segmenter,
):
    unsegmented_text = get_text(
        unsegmented,
        json_unsegmented_transcript_key,
    )
    segmented_text = get_text(
        segmented,
        json_segmented_transcript_key,
    )
    alignment = pairwise2.align.globalcx(
        segmented_text,
        unsegmented_text,
        match_fn_mwer if mwer_segmenter else match_fn,
        one_alignment_only=True,
        gap_char=[None],
    )[0]
    segmented_text = get_text(
        segmented,
        json_segmented_transcript_key,
        keep_lines=True,
    )

    segmented_new = []
    last = 0
    for original_line in segmented_text:
        s, ns = [], []
        while len(s) < len(original_line):
            if alignment.seqA[last]:
                s.append(alignment.seqA[last])
            if alignment.seqB[last]:
                ns.append(alignment.seqB[last])
            last += 1
        segmented_new.append(
            {
                json_unsegmented_transcript_key: " ".join(ns),
                "prediction_length": len(ns),
            }
        )
    return segmented_new


def read_json(file):
    lines = []
    with open(file) as file:
        for line in file:
            line = json.loads(line)
            lines.append(line)
    return lines


def join_metadata(unsegmented, keys):
    metadata = {}
    for key in keys:
        data = []
        for line in unsegmented:
            data += line[key]
        metadata[key] = data
    return metadata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--unsegmented",
        type=str,
        default="unsegmented.log",
    )
    parser.add_argument(
        "--segmented",
        type=str,
        default="segmented.log",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="output.log",
    )
    parser.add_argument(
        "--json-unsegmented-transcript-key", type=str, default="prediction"
    )
    parser.add_argument(
        "--json-segmented-transcript-key", type=str, default="reference"
    )
    parser.add_argument(
        "--json-unsegmented-metadata-keys", nargs="+", default=["delays", "elapsed"]
    )
    parser.add_argument(
        "--json-segmented-metadata-keys",
        nargs="+",
        default=["reference", "source_length", "reference_length"],
    )
    parser.add_argument("--mwer-segmenter", action="store_true")

    args = parser.parse_args()

    segmented = read_json(args.segmented)
    unsegmented = read_json(args.unsegmented)
    new_segmented = segment(
        unsegmented,
        segmented,
        args.json_unsegmented_transcript_key,
        args.json_segmented_transcript_key,
        args.mwer_segmenter,
    )
    unsegmented_metadata = join_metadata(
        unsegmented, args.json_unsegmented_metadata_keys
    )

    total = 0
    for s, ns in zip(segmented, new_segmented):
        nsl = ns["prediction_length"]
        for key in args.json_unsegmented_metadata_keys:
            ns[key] = unsegmented_metadata[key][total : total + nsl]
        total += nsl

        for key in args.json_segmented_metadata_keys:
            ns[key] = s[key]

    with open(args.output, "w") as output:
        for line in new_segmented:
            output.write(json.dumps(line, ensure_ascii=False))
            output.write("\n")


if __name__ == "__main__":
    main()
