"""Meaningful control-flow and honest-transcript regression checks; untimed."""
import copy
import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import protocol as P
from session import SequentialSession, SessionError

H = ((1, 0, 1, 0), (0, 1, 0, 1))
E = (1, 1, 0, 0)


def context():
    return P.Context(5, 4, 2, (2, 2, 0, 0, 0), H, (1, 1))


class StateRuleTests(unittest.TestCase):
    def test_each_invalid_challenge_consumes_state(self):
        for bad in (3, -1, True, None, '0'):
            with self.subTest(challenge=bad):
                A, state = P.Prover(context(), E).commit()
                with self.assertRaises(ValueError):
                    P.respond(state, bad)
                self.assertTrue(state.used)
                with self.assertRaises(ValueError):
                    P.respond(state, 0)

    def test_three_valid_challenges_and_reuse_rejection(self):
        for b in (0, 1, 2):
            ctx = context()
            A, state = P.Prover(ctx, E).commit()
            response = P.respond(state, b)
            self.assertTrue(P.verify(ctx, A, b, response))
            for again in (0, 1, 2, 3):
                with self.assertRaises(ValueError):
                    P.respond(state, again)

    def test_same_honest_transcripts_as_historical_implementation(self):
        # Use the original workspace or the bundled frozen source snapshot.
        path = Path(__file__).resolve().parents[1] / 'performance/length_sweep/implementation/protocol.py'
        if not path.is_file():
            path = Path(__file__).resolve().parents[1] / 'historical_evidence/frozen_implementation/protocol.py'
        if not path.is_file():
            self.skipTest('no historical source snapshot available')
        spec = importlib.util.spec_from_file_location('historical_protocol', path)
        old = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = old
        spec.loader.exec_module(old)
        ctx = context()
        A, state = P.Prover(ctx, E).commit()
        old_ctx = old.Context(5, 4, 2, (2, 2, 0, 0, 0), H, (1, 1))
        for b in (0, 1, 2):
            legacy_state = old.State(state.sigma, state.u, state.z, state.e, state.salts)
            new_response = P.respond(copy.deepcopy(state), b)
            old_response = old.respond(legacy_state, b)
            self.assertEqual(new_response, old_response)
            self.assertTrue(old.verify(old_ctx, A, b, new_response))

    def test_full_round_count_required(self):
        session = SequentialSession(P.Prover(context(), E), context(), 3)
        for b in (0, 1, 2):
            self.assertFalse(session.accepted)
            session.start_round()
            session.answer(b)
            self.assertTrue(session.finish_round())
        self.assertTrue(session.accepted)
        self.assertEqual(session.rounds_done, 3)
        with self.assertRaises(SessionError):
            session.start_round()
        self.assertTrue(session.accepted)

    def test_invalid_challenge_closes_whole_session(self):
        session = SequentialSession(P.Prover(context(), E), context(), 2)
        session.start_round()
        state = session._state
        with self.assertRaises(ValueError):
            session.answer(3)
        self.assertEqual(session.phase, 'closed')
        self.assertTrue(state.used)
        with self.assertRaises(SessionError):
            session.answer(0)

    def test_abort_and_out_of_order_invalidate_pending_state(self):
        for action in ('abort', 'duplicate_commit', 'premature_verify'):
            with self.subTest(action=action):
                session = SequentialSession(P.Prover(context(), E), context(), 2)
                session.start_round()
                state = session._state
                if action == 'abort':
                    session.abort('transport disconnected')
                else:
                    with self.assertRaises(SessionError):
                        session.start_round() if action == 'duplicate_commit' else session.finish_round()
                self.assertTrue(state.used)
                self.assertEqual(session.phase, 'closed')

    def test_duplicate_answer_cannot_expose_second_branch(self):
        session = SequentialSession(P.Prover(context(), E), context(), 2)
        session.start_round()
        session.answer(0)
        with self.assertRaises(SessionError):
            session.answer(1)
        self.assertEqual(session.phase, 'closed')

    def test_rejected_response_terminates_authentication(self):
        session = SequentialSession(P.Prover(context(), E), context(), 2)
        session.start_round()
        response = session.answer(0)
        response['salts'] = (b'X' * 32, b'Y' * 32)
        self.assertFalse(session.finish_round(response))
        self.assertFalse(session.accepted)
        self.assertEqual(session.phase, 'closed')

    def test_execution_errors_terminate(self):
        session = SequentialSession(P.Prover(context(), E), context(), 2)
        with patch.object(session.prover, 'commit', side_effect=RuntimeError('random source unavailable')):
            with self.assertRaises(RuntimeError):
                session.start_round()
        self.assertEqual(session.phase, 'closed')

    def test_round_count_cannot_be_negotiated_as_unvalidated_input(self):
        for bad in (0, -1, True, '219'):
            with self.assertRaises(ValueError):
                SequentialSession(P.Prover(context(), E), context(), bad)


if __name__ == '__main__':
    unittest.main(verbosity=2)
