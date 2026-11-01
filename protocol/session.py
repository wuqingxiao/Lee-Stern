"""Local sequential session control; no network framing or timeout transport."""
import secrets
import protocol as P


class SessionError(ValueError):
    pass


class SequentialSession:
    """Own one pending commitment and terminate on invalid order/error/rejection.

    This local harness models both endpoints. It is not a network authenticator.
    The caller supplies a trusted round count and handles any transport abort by
    calling abort(). Honest challenges can be generated with random_challenge().
    """

    def __init__(self, prover, verifier_context, rounds):
        if type(rounds) is not int or rounds < 1:
            raise ValueError('trusted positive round count required')
        self.prover = prover
        self.ctx = verifier_context
        self.rounds = rounds
        self.rounds_done = 0
        self._phase = 'commit'
        self._state = self._A = self._b = self._response = None
        self.reason = None

    @property
    def phase(self):
        return self._phase

    @property
    def accepted(self):
        return self._phase == 'complete'

    def abort(self, reason='session aborted'):
        if self._phase == 'complete':
            return
        if self._state is not None:
            self._state.used = True
        self._phase = 'closed'
        self.reason = reason
        self._state = self._A = self._b = self._response = None

    def _require(self, phase):
        if self._phase in ('closed', 'complete'):
            raise SessionError('session is terminal')
        if self._phase != phase:
            self.abort('invalid message order')
            raise SessionError('invalid message order')

    def start_round(self):
        self._require('commit')
        try:
            self._A, self._state = self.prover.commit()
        except Exception:
            self.abort('commitment generation failed')
            raise
        self._phase = 'answer'
        return self._A

    def answer(self, challenge):
        self._require('answer')
        try:
            response = P.respond(self._state, challenge)
        except Exception:
            self.abort('invalid challenge or response error')
            raise
        self._b, self._response = challenge, response
        self._phase = 'verify'
        return response

    def random_challenge(self):
        self._require('answer')
        try:
            return secrets.randbelow(3)
        except Exception:
            self.abort('challenge generation failed')
            raise

    def finish_round(self, response=None):
        self._require('verify')
        supplied = self._response if response is None else response
        try:
            ok = P.verify(self.ctx, self._A, self._b, supplied)
        except Exception:
            self.abort('verification error')
            raise
        if not ok:
            self.abort('verification rejected')
            return False
        self.rounds_done += 1
        self._state = self._A = self._b = self._response = None
        self._phase = 'complete' if self.rounds_done == self.rounds else 'commit'
        return True
