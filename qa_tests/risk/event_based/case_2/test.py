# Copyright (c) 2010-2012, GEM Foundation.
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

import os

from nose.plugins.attrib import attr as noseattr

from qa_tests import risk
from tests.utils import helpers

from openquake.engine.db import models

# FIXME(lp). This is a regression test. Data has not been validated
# by an alternative reliable implemantation


class EventBasedRiskCase2TestCase(risk.BaseRiskQATestCase):
    cfg = os.path.join(os.path.dirname(__file__), 'job.ini')

    @noseattr('qa', 'risk', 'event_based')
    def test(self):
        self._run_test()

    def hazard_id(self):
        job = helpers.get_hazard_job(
            helpers.get_data_path("event_based_hazard/job.ini"))
        gmf_coll = helpers.create_gmf_from_csv(job, os.path.join(
            os.path.dirname(__file__), 'gmf.csv'))

        return gmf_coll.output.id

    def actual_data(self, job):
        data = ([curve.poes
                for curve in models.LossCurveData.objects.filter(
                    loss_curve__output__oq_job=job,
                    loss_curve__aggregate=False,
                    loss_curve__insured=False).order_by('asset_ref')] +
                [curve.loss_ratios
                for curve in models.LossCurveData.objects.filter(
                    loss_curve__output__oq_job=job,
                    loss_curve__aggregate=False,
                    loss_curve__insured=False).order_by('asset_ref')] +
                [curve.losses
                for curve in models.AggregateLossCurveData.objects.filter(
                    loss_curve__output__oq_job=job,
                    loss_curve__aggregate=True,
                    loss_curve__insured=False)] +
                [[el.aggregate_loss
                 for el in models.EventLoss.objects.filter(
                output__oq_job=job).order_by('-aggregate_loss')[0:10]]])
        return data

    def expected_data(self):

        poes_1 = [1.0, 1.0, 0.981684361111, 0.864664716763, 0.864664716763,
                  0.632120558829, 0.632120558829, 0.632120558829, 0.0]

        poes_2 = [1.0, 1.0, 1.0, 0.999999999986, 0.999997739671,
                  0.999876590196, 0.999088118034, 0.981684361111, 0.0]

        poes_3 = [1.0, 1.0, 1.0, 0.99995460007, 0.997521247823,
                  0.981684361111, 0.864664716763, 0.864664716763, 0.0]

        losses_1 = [0.0, 0.0112395242406, 0.0224790484812, 0.0337185727218,
                    0.0449580969624, 0.056197621203, 0.0674371454436,
                    0.0786766696842, 0.0899161939248]

        losses_2 = [0.0, 0.0019748458149, 0.0039496916298, 0.0059245374447,
                    0.00789938325961, 0.00987422907451, 0.0118490748894,
                    0.0138239207043, 0.0157987665192]

        losses_3 = [0.0, 0.00605092131392, 0.0121018426278, 0.0181527639418,
                    0.0242036852557, 0.0302546065696, 0.0363055278835,
                    0.0423564491975, 0.0484073705114]

        expected_aggregate_losses = [0.0, 41.9450370893, 83.8900741787,
                                     125.835111268, 167.780148357,
                                     209.725185447, 251.670222536,
                                     293.615259625, 335.560296715]

        expected_event_loss_table = [335.560296714705, 219.63382977226,
                                     161.101434197127, 117.587180097578,
                                     116.488051808864, 107.949974628949,
                                     106.171317215111, 105.729211465077,
                                     98.34853475773, 93.4503401489587]

        return [poes_1, poes_2, poes_3, losses_1, losses_2, losses_3,
                expected_aggregate_losses, expected_event_loss_table]

    def actual_xml_outputs(self, job):
        """
        do not check file outputs
        """
        return []

    def expected_outputs(self):
        return []
