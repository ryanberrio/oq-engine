# Copyright (c) 2013, GEM Foundation.
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake.  If not, see <http://www.gnu.org/licenses/>.


from nose.plugins.attrib import attr as noseattr
from qa_tests import risk


class EventBaseQATestCase(risk.CompleteTestCase, risk.FixtureBasedQATestCase):
    hazard_calculation_fixture = ("QA (regression) test for Risk Event "
                                  "Based from Stochastic Event Set")

    @noseattr('qa', 'risk', 'event_based')
    def test(self):

        expected_losses = [
            124.194084993, 10.551647246, 9.28007570766, 5.24747652954,
            2.65924996464, 1.95625047528, 1.04670143538, 0.917858366911,
            0.0795391520762, 0.016769321535, 0.0116411543827]

        losses = self._run_test().output_set.get(
            output_type="event_loss").event_loss
        # uncomment the following line if the test fails and you want to see
        # the numbers (this will happen frequently)
        # print [l.aggregate_loss for l in losses]
        for event_loss, expected in zip(losses, expected_losses):

            self.assertAlmostEqual(
                expected, event_loss.aggregate_loss,
                msg="loss for rupture %r is %s (expected %s)" % (
                    event_loss.rupture.tag, event_loss.aggregate_loss,
                    expected))
