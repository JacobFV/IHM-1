# Canonical execution baseline

The approved implementation starts at `f2dbcefe5e3ebd4ad2bb053d6fb36597daf07a45` on `feat/integrated-human`. The accompanying machine-readable baseline is `data/derived/audits/execution-baseline-20260905.json`; the repeatable audit is `scripts/audit_body_coverage.py`.

At this revision the canonical trajectory has 301 recorded frames and seven moving entity transforms: five lung lobes and two cardiac cavity surfaces. The anatomy has 2,408 representations. These are different denominators and neither establishes whole-body dynamic coverage.

The committed respiration and peripheral solvers are retained as reduced prototypes to integrate and test. The fixed workspace is retained for browser verification. The peripheral prototype's implicit touch-gated cortical motor readout will be removed from default body control: sensory arrival does not justify arbitrary muscle recruitment. The coarse affine mechanics remains explicitly bounded; it cannot qualify as reference volumetric contact mechanics.

The running loopback server initially used PID 1448779. Its process started before these changes; final deployment must restart it after checking for active jobs. Existing generated source bytes and old run receipts are preserved. A current source manifest does not retroactively certify an older trajectory.

Scientific gaps remain recorded independently from implementation and numerical checks. In particular source geometry does not supply microvascular connectivity, material interfaces, innervation, individualized coefficients, or a validated regenerative response law.
