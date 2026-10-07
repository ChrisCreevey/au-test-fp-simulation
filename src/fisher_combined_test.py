#!/usr/bin/env python3
"""
Compare the real delta-logL against each domain's simulated null
distribution, and combine the resulting per-domain p-values into one
dataset-wide significance test.

Per domain: a one-sided empirical/Monte Carlo permutation p-value --
P(null delta-logL >= real delta-logL), with the standard +1/+1 continuity
correction (Davison & Hinkley) so the p-value is never exactly 0 even when
the real value is more extreme than every null replicate:

    p = (n_null_as_extreme + 1) / (n_valid_replicates + 1)

This tests, per domain: is the real data's evidence against monophyly
bigger than pure simulation noise under the true null model would
produce? Valid at any replicate count, but coarse-grained at N=10 (only
11 possible values, k/11 for k=1..11) -- fine for aggregate use, not
precise enough to call an individual domain significant on its own.

Combined across domains: Fisher's method,
    X = -2 * sum(ln(p_i))  ~  chi-squared with 2k degrees of freedom
(k = number of domains combined), testing whether the real data is
*systematically* more extreme than the null simulation predicts across
the whole dataset -- the aggregate answer to "is 41% surprising?", not
just a count of individually-significant domains. Main assumption:
independence across domains (reasonable as an approximation across
different gene families, not exact).

The chi-squared survival function for EVEN degrees of freedom (always
true here, df=2k) has an exact closed form (no scipy needed):
    P(X > x) = exp(-x/2) * sum_{i=0}^{k-1} (x/2)^i / i!
computed in log-space below, since k can be in the thousands and the
naive form overflows/underflows badly at that scale.

Usage:
  python3 src/fisher_combined_test.py --report data/reports/failure_report_TIMESTAMP.tsv \
      [--out-tsv data/reports/fisher_test_detail.tsv] [--email you@example.com]
"""
import argparse
import csv
import math
import re
import smtplib
import sys
from collections import defaultdict
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path


def permutation_p_upper(real_delta, null_deltas):
    """One-sided permutation p-value: P(null >= real), +1/+1 corrected."""
    n_as_extreme = sum(1 for v in null_deltas if v >= real_delta)
    return (n_as_extreme + 1) / (len(null_deltas) + 1)


def log_chi2_sf_even_df(x, k):
    """log(P(chi2_{2k} > x)), computed via the exact Poisson-CDF relation
    for even degrees of freedom, in log-space to avoid overflow/underflow
    when k (domain count) or x (Fisher's statistic) is large.
    P(chi2_2k > x) = exp(-x/2) * sum_{i=0}^{k-1} (x/2)^i / i!
    """
    if k <= 0:
        return 0.0  # no domains combined -> p = 1 (log 1 = 0)
    if x <= 0:
        return 0.0  # statistic is 0 -> p = 1

    half_x = x / 2.0
    log_half_x = math.log(half_x)

    # log of each term i*log(x/2) - lgamma(i+1), then logsumexp over i=0..k-1
    log_terms = [i * log_half_x - math.lgamma(i + 1) for i in range(k)]
    m = max(log_terms)
    log_sum = m + math.log(sum(math.exp(t - m) for t in log_terms))

    return -half_x + log_sum


def fishers_combined(pvalues):
    """Returns (statistic, log10_pvalue, k). log10_pvalue is reported
    (not the raw p) since p can be far smaller than floating-point can
    represent once hundreds/thousands of domains are combined."""
    k = len(pvalues)
    if k == 0:
        return None, None, 0
    stat = -2.0 * sum(math.log(p) for p in pvalues)
    log_p = log_chi2_sf_even_df(stat, k) / math.log(10)
    return stat, log_p, k


def collect_domain_pvalues(report_path):
    """Returns:
      by_domain: {domain_slot: [(family, p_value), ...]} -- for Fisher's method
      by_family: {family: {domain_slot: (real_delta, n_null, p_value)}} --
        per-family detail, each domain in its own slot, matching the wide
        (one row per family) convention every other report in this
        pipeline uses.
    """
    by_domain = defaultdict(list)
    by_family = defaultdict(dict)

    with open(report_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        domain_slots = sorted({
            m.group(1) for col in reader.fieldnames
            if (m := re.match(r'^Domain(\d+)_Real_DeltaLogL$', col))
        }, key=int)

        rep_cols = {}
        for d in domain_slots:
            rep_cols[d] = sorted(
                (int(m.group(1)), col) for col in reader.fieldnames
                if (m := re.match(rf'^Domain{d}_DeltaLogL_Rep(\d+)$', col))
            )

        for row in reader:
            fam = row["Family_ID"]
            for d in domain_slots:
                real_raw = row.get(f"Domain{d}_Real_DeltaLogL", "NA")
                if real_raw == "NA" or real_raw == "":
                    continue
                try:
                    real_delta = float(real_raw)
                except ValueError:
                    continue

                null_deltas = []
                for _, col in rep_cols[d]:
                    v = row.get(col, "NA")
                    if v != "NA" and v != "":
                        try:
                            null_deltas.append(float(v))
                        except ValueError:
                            pass

                if not null_deltas:
                    continue

                p = permutation_p_upper(real_delta, null_deltas)
                by_domain[d].append((fam, p))
                by_family[fam][d] = (real_delta, len(null_deltas), p)

    return by_domain, by_family


def send_email(to_addr, subject, body, attachment_path=None):
    msg = MIMEMultipart()
    msg["Subject"] = subject
    msg["From"] = "au-test-sim@localhost"
    msg["To"] = to_addr
    msg.attach(MIMEText(body, "plain"))
    if attachment_path:
        with open(attachment_path, "rb") as f:
            part = MIMEApplication(f.read(), Name=Path(attachment_path).name)
        part["Content-Disposition"] = f'attachment; filename="{Path(attachment_path).name}"'
        msg.attach(part)
    with smtplib.SMTP("localhost") as smtp:
        smtp.send_message(msg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True, help="generate_failure_report.py output TSV")
    ap.add_argument("--out-tsv", default=None, help="Optional: write per-domain p-values to this file")
    ap.add_argument("--email", default=None, help="If given, email the summary to this address")
    args = ap.parse_args()

    by_domain, by_family = collect_domain_pvalues(args.report)

    if not by_domain:
        print("No domains with both a real delta-logL and null replicate values found.", file=sys.stderr)
        sys.exit(1)

    lines = [
        "Fisher's combined test: real delta-logL vs. simulated null distribution",
        "",
        "Per-domain test: one-sided permutation p-value, P(null delta-logL >= real),",
        "+1/+1 continuity corrected. Combined: Fisher's method",
        "(-2 * sum(ln(p_i)) ~ chi-squared, df=2k) across all domains with both a",
        "real delta-logL and null replicates available.",
        "",
    ]

    pooled = []
    for d in sorted(by_domain.keys(), key=int):
        pvals = [p for _, p in by_domain[d]]
        pooled.extend(pvals)
        stat, log10_p, k = fishers_combined(pvals)
        n_sig_05 = sum(1 for p in pvals if p <= 0.05)
        lines.append(f"Domain{d}: {k} domains tested")
        lines.append(f"  Individually significant at p<=0.05: {n_sig_05}/{k} ({100*n_sig_05/k:.1f}%)")
        lines.append(f"  Fisher's statistic: {stat:.1f} (df={2*k})")
        lines.append(f"  Combined log10(p-value): {log10_p:.1f}  (p {'< 1e-300' if log10_p < -300 else f'~ 10^{log10_p:.1f}'})")
        lines.append("")

    stat, log10_p, k = fishers_combined(pooled)
    n_sig_05 = sum(1 for p in pooled if p <= 0.05)
    lines.append(f"POOLED across all domain slots: {k} domain-tests")
    lines.append(f"  Individually significant at p<=0.05: {n_sig_05}/{k} ({100*n_sig_05/k:.1f}%)")
    lines.append(f"  Fisher's statistic: {stat:.1f} (df={2*k})")
    lines.append(f"  Combined log10(p-value): {log10_p:.1f}  (p {'< 1e-300' if log10_p < -300 else f'~ 10^{log10_p:.1f}'})")

    body = "\n".join(lines)
    print(body)

    if args.out_tsv:
        domain_slots = sorted(by_domain.keys(), key=int)
        header = ["Family_ID"]
        for d in domain_slots:
            header += [f"Domain{d}_Real_DeltaLogL", f"Domain{d}_N_Null_Replicates",
                       f"Domain{d}_Permutation_P_Upper"]

        with open(args.out_tsv, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(header)
            for fam in sorted(by_family.keys(), key=lambda s: int(s[1:]) if s[1:].isdigit() else 0):
                row = [fam]
                for d in domain_slots:
                    if d in by_family[fam]:
                        real_delta, n_null, p = by_family[fam][d]
                        row += [real_delta, n_null, p]
                    else:
                        row += ["NA", "NA", "NA"]
                writer.writerow(row)
        print(f"\nPer-family, per-domain detail written to: {args.out_tsv}")

    if args.email:
        send_email(args.email, "Fisher's combined test: real vs. null delta-logL", body, args.out_tsv)
        print(f"\nEmailed to {args.email}", file=sys.stderr)


if __name__ == "__main__":
    main()
