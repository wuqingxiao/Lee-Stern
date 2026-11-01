# Known issues in the archived snapshot

This file prevents the research archive from being mistaken for a final security certification.

## PKP Wagner partition

`experiments/p2_pkp/pkp_price_v2.py` currently uses four blocks with `floor((k+u)/4)` coordinates per block. For the main instance `k+u=779`, this covers 776 of 779 coordinates and leaves three coordinates outside the partition. The recorded `2^711.5` value therefore cannot be presented as the cost of a complete disjoint-block Wagner attack. Treat it as an incomplete or simplified model value until the partition is corrected and the result is regenerated.

## Cost-unit comparison

KMP bit-operation estimates, Wagner list-operation exponents, BKZ core-SVP estimates, and one-call ISD/FCSD figures are different accounting units. They must not be sorted as one universal “cheapest attack” table without a justified conversion.

## Delivery metadata

Some working-directory delivery metadata was generated before the final audit and may contain a stronger approval flag than `security/FINAL_PARAMETER_DECISION.*`. The final decision files and `audit/STATUS_AND_LIMITATIONS.md` are authoritative for this snapshot.
