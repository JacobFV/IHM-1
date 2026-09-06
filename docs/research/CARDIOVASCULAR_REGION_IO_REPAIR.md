# Cardiovascular region metadata repair

This isolated repair preserves an existing Circuit schema field. It does not recover physiological history or promote a production library. The preceding signed cardiovascular experiment is documented in `SIGNED_CARDIOVASCULAR_READER_AND_PRIOR.md`; its missing regional resistance response remains a retained failure until an independently verified repaired experiment succeeds.

The pinned `src/io/cdm/Circuit.cpp` omits `CardiovascularRegion` in both fluid-path serialization directions, although `Circuit.xsd` and generated `Circuit.hxx` already provide the optional field. `Cardiovascular::SetUp` uses that field to populate regional resistance caches. An unlabeled legacy snapshot therefore produces empty muscle resistance caches even when the signed demand reader receives its input correctly.

`patch_cardiovascular_region_io.py` adds explicit mappings for all five native/schema regions. Both directions clear a reused destination field before restoring an optional value. Unexpected enum values reject. The builder replaces only `io/cdm/Circuit.cpp.o`, verifies inherited object identities, and links an isolated library. Original pinned headers retain one include identity; a first attempted build with copied pragma-once headers failed because `SECircuit.inl` also includes their original paths. That failed directory is retained as `whole_body_integrity_cardio_region_io_v1`.

Legacy migration is a separate explicit operation. It accepts only the exact retained exertion snapshot hash and pinned controller source. Parsing native construction order yields 41 assignments, six subsequent deletions, and 35 final region labels. Four deleted paths are absent in the retained graph; two kidney inlet path names are recreated without region labels and deliberately remain unlabeled. Every final assigned path must exist. The migration records the 348-path graph identity and inserts only the 35 metadata tags, in existing schema order after resistance baseline. Removing the inserted tags must reproduce every original byte. Missing paths, changed source, an already labeled snapshot, or another state identity reject.

This is deterministic topology reconstruction. It changes which resistance paths participate in existing regional laws after loading; it is therefore not a claim that unlabeled and repaired trajectories should match. Subsequent zero-demand parity must compare the experimental prior against the same region-repaired baseline. No dynamic scalar, private physiological latch, or elapsed history is reconstructed.

Validation entry points:

- `scripts/verify_cardiovascular_region_io.py`: pinned source, final topology, exact nonmetadata preservation, and receipt guards.
- `scripts/verify_cardiovascular_region_io_native.py --run --variant NAME`: isolated native five-region DTO roundtrip, retained resistance baselines, and clearing of absent values in reused native/serialized destinations. This fixture does not load or advance a patient.

The nervous sleep variant has a different paired CDM/schema inventory. This builder explicitly rejects that branch until a composition recipe verifies and retains its complete paired object/library identities; it must not replace that branch with an older CDM library.

The corrected isolated library `whole_body_integrity_cardio_region_io_v2` built in 2.248 seconds with maximum child RSS 360076 KiB (library SHA-256 `1d412022ed9ae1be213098cada40d9975b85e2ff20649d17bda0dc39ef4fd80e`). Its one replaced object SHA-256 is `06500cdb6152b3f0651e6b42d7e77f1aa8f86622c506f19c4fc9f50e939c3d37`.

The first standalone roundtrip link demonstrated that `io::Circuit` methods are private library symbols. The fixture therefore links the exact corrected Circuit object and inherited Property object directly, with their manifest hashes verified, alongside the normal public library dependencies. Symbol inspection identified precisely two hidden undefined references from Circuit, both Property enum conversions; Property has no further hidden undefined references. This does not alter library symbol visibility and is explicitly an object-level IO test, distinct from a patient load using the linked shared library. Failed compile/link receipts are retained.

The object-level native roundtrip passed in 1.138 seconds with maximum child RSS 285236 KiB, retained at `data/derived/audits/cardiovascular-region-io-native-o_oglcjn`. All five regions, resistance baselines, and absent-field clearing passed. The inherited Property object also requires an explicit Xerces library link. No patient advance occurred. Actual resistance feedback after metadata restoration still requires the separately queued shared-library experiment.
