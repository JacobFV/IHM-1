"""Independent copies of internal snapshot data, without Python tree traversal."""
import pickle


def clone_snapshot_data(value):
    """Clone an in-memory, pickle-compatible internal data graph.

    Preserve types, repeated references, and cycles while isolating mutable
    snapshot data from its owner. This is not a serialized-input interface:
    the only bytes decoded are produced immediately here from ``value``.
    Never replace this round trip with loading caller-supplied bytes or files.
    Native engine handles and arbitrary external objects are not snapshot data.
    """
    return pickle.loads(pickle.dumps(value, protocol=5))
