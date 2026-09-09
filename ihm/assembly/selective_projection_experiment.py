"""Experimental internal subset projection; production runtime does not import this.

Public advance, state and snapshot retain complete plain dictionaries. Only the
explicit advance_observation method returns a subset. Native force routing and
integration are unchanged. A private captured endpoint permits a subsequent
full snapshot without observing a newer native state.
"""
import copy

from .articulated import ArticulatedBodyPlant
from .snapshot_data import clone_snapshot_data


class SelectiveProjectionPlant(ArticulatedBodyPlant):
    @property
    def state(self):
        if self._full_frame is None:
            self._full_frame = self._project(self._projection_endpoint, self._projection_work)
        return self._full_frame

    @state.setter
    def state(self, frame):
        self._full_frame = frame
        self._projection_endpoint = None
        self._projection_work = None

    def _project_observation(self, native, work, entity_ids):
        # A normal private registration copy has a smaller specs dictionary;
        # reuse the exact production projection arithmetic and frame assembly.
        # The actual registration, its force routing and its manifest stay full.
        registration = copy.copy(self.registration)
        registration.specs = {key: value for key, value in self.registration.specs.items()
                              if key in entity_ids}
        view = object.__new__(ArticulatedBodyPlant)
        view.registration = registration
        view.garments = self.garments
        return ArticulatedBodyPlant._project(view, native, work)

    def advance_observation(self, dt_s, forces=(), actuation=None, *, entity_ids):
        ids = frozenset(entity_ids)
        if ids - self.registration.specs.keys():
            raise ValueError('Unknown selective canonical entity')
        old = self.native.snapshot()
        mapped = []
        for force in forces:
            if set(force) != {'id', 'point_m', 'force_n'}:
                raise ValueError('Canonical force requires id, point_m, force_n')
            mapped.append(self.registration.force(force['id'], force['point_m'], force['force_n'], old))
        checkpoint = self.native.checkpoint()
        garment_checkpoint = None if self.garments is None else self.garments.checkpoint()
        cached = self._full_frame, self._projection_endpoint, self._projection_work
        try:
            result = (self.native.advance(dt_s, mapped, actuation) if self.garments is None else
                      self.garments.advance(self.native, dt_s, mapped, actuation)[0])
            work = result['positive_active_fiber_work_j'] - old['positive_active_fiber_work_j']
            observation = self._project_observation(result, work, ids)
            # No native-result aliases escape or become future snapshot inputs.
            endpoint = clone_snapshot_data(result)
            answer = clone_snapshot_data(observation)
            self._full_frame = None
            self._projection_endpoint = endpoint
            self._projection_work = work
            return answer
        except BaseException:
            self.native.restore(checkpoint)
            if self.garments is not None:
                self.garments.restore(garment_checkpoint)
            self._full_frame, self._projection_endpoint, self._projection_work = cached
            raise
        finally:
            self.native.release(checkpoint)


class SubstepObservationOwner:
    """Explicit experimental internal interface; never a public snapshot proxy.

    Existing exchange helpers only call advance. Complete their final subset
    with owner.snapshot() within the enclosing exchange rollback boundary,
    retaining separately aggregated work/audits. A production integration should
    instead use normal owner.advance on the final substep so full projection is
    inside that substep's transaction. Production does not use this adapter.
"""
    def __init__(self, owner, entity_ids):
        self.owner = owner
        self.entity_ids = frozenset(entity_ids)

    def advance(self, dt_s, forces=(), actuation=None):
        return self.owner.advance_observation(dt_s, forces, actuation, entity_ids=self.entity_ids)
