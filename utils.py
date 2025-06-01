# utils.py
import io
from Bio import SeqIO
from Bio import pairwise2

# --- Utility Functions ---

def load_sequence(handle_or_path):
    """
    Read the first sequence from a file path or uploaded FASTA (file-like).
    Returns the sequence string or None on failure.
    """
    try:
        if isinstance(handle_or_path, str):
            with open(handle_or_path, 'r') as handle:
                for rec in SeqIO.parse(handle, 'fasta'):
                    return str(rec.seq)
            return None
        else:
            handle_or_path.seek(0)
            data = handle_or_path.read()
            if isinstance(data, bytes):
                data = data.decode('utf-8')
            text_handle = io.StringIO(data)
            for rec in SeqIO.parse(text_handle, 'fasta'):
                return str(rec.seq)
            return None
    except Exception as e:
        raise RuntimeError(f"Failed to load sequence: {e}")


def extract_s_gene(full_seq, start=21563, end=25384):
    """
    Extract the Spike (S) gene slice from a full-genome sequence (1-based inclusive).
    Returns the S-gene sequence string.
    """
    slice_start = start - 1
    slice_end = end
    if not full_seq or len(full_seq) < slice_end:
        raise ValueError("Sequence too short for S-gene extraction.")
    return full_seq[slice_start:slice_end]


def align_sequences(seq1, seq2,
                    match_score=2, mismatch_score=-1,
                    gap_open_penalty=-0.5, gap_extend_penalty=-0.1):
    """
    Perform global alignment between two sequences using pairwise2.globalms.
    Returns the best alignment tuple:
        (aligned_seq1, aligned_seq2, score, begin, end)
    """
    alignments = pairwise2.align.globalms(
        seq1, seq2,
        match_score, mismatch_score,
        gap_open_penalty, gap_extend_penalty
    )
    return alignments[0] if alignments else None


def find_mutations(reference_aligned_seq, variant_aligned_seq, original_reference_seq):
    """
    Identify mutations (substitutions, insertions, deletions)
    from two aligned sequences. Returns list of mutation dicts:
      - Substitution: {position_ref, type, ref_base, var_base}
      - Insertion:   {position_ref_after, type, inserted_bases}
      - Deletion:    {position_ref, type, deleted_ref_base}
    Positioning uses 1-based original reference coordinates.
    """
    mutations_list = []
    original_ref_pos = 0
    for i in range(len(reference_aligned_seq)):
        ref_base = reference_aligned_seq[i]
        var_base = variant_aligned_seq[i]
        if ref_base != '-':
            original_ref_base = original_reference_seq[original_ref_pos] if original_ref_pos < len(original_reference_seq) else None
        else:
            original_ref_base = None

        if ref_base != var_base:
            # substitution
            if ref_base != '-' and var_base != '-':
                mutations_list.append({
                    "position_ref": original_ref_pos + 1,
                    "type": "Substitution",
                    "ref_base": original_ref_base,
                    "var_base": var_base
                })
            # insertion
            elif ref_base == '-' and var_base != '-':
                mutations_list.append({
                    "position_ref_after": original_ref_pos,
                    "type": "Insertion",
                    "inserted_bases": var_base
                })
            # deletion
            elif ref_base != '-' and var_base == '-':
                mutations_list.append({
                    "position_ref": original_ref_pos + 1,
                    "type": "Deletion",
                    "deleted_ref_base": original_ref_base
                })
        if ref_base != '-':
            original_ref_pos += 1
    return mutations_list


def count_matching_mutations(mutations_list_patient, mutations_list_variant_profile):
    """
    Count how many mutation events from patient list appear in variant profile list.
    Uses tuple(sorted(dict.items())) as signature for comparison.
    """
    def sig(m):
        return tuple(sorted(m.items()))
    patient_sigs = {sig(m) for m in mutations_list_patient}
    variant_sigs = {sig(m) for m in mutations_list_variant_profile}
    return len(patient_sigs & variant_sigs)


def classify_variant(patient_mutations, variant_profiles, wuhan_like_threshold=3):
    """
    Given patient mutation list and dict of variant->mutation_list,
    return (identified_variant, match_count).
    Applies threshold: if no patient mutations -> 'Wuhan-like';
    if max matches < threshold -> 'Wuhan-like/Unknown'; else select best or tie.
    """
    if not patient_mutations:
        return ("Wuhan-like (No mutations)", 0)
    scores = {var: count_matching_mutations(patient_mutations, profile)
              for var, profile in variant_profiles.items()}
    best_var = max(scores, key=scores.get)
    best_score = scores[best_var]
    if best_score < wuhan_like_threshold:
        return ("Wuhan-like or Other (Low match)", best_score)
    ties = [var for var, sc in scores.items() if sc == best_score]
    if len(ties) > 1:
        return (" or ".join(ties) + " (tie)", best_score)
    return (best_var, best_score)