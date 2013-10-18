# Copyright (c) 2010-2013, GEM Foundation.
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

"""
Core calculator functionality for computing stochastic event sets and ground
motion fields using the 'event-based' method.

Stochastic events sets (which can be thought of as collections of ruptures) are
computed iven a set of seismic sources and investigation time span (in years).

For more information on computing stochastic event sets, see
:mod:`openquake.hazardlib.calc.stochastic`.

One can optionally compute a ground motion field (GMF) given a rupture, a site
collection (which is a collection of geographical points with associated soil
parameters), and a ground shaking intensity model (GSIM).

For more information on computing ground motion fields, see
:mod:`openquake.hazardlib.calc.gmf`.
"""

import copy
import random

import openquake.hazardlib.imt
import numpy.random

from django.db import transaction
from openquake.hazardlib.calc import filters
from openquake.hazardlib.calc import gmf
from openquake.hazardlib.calc import stochastic

from openquake.engine import writer
from openquake.engine.calculators.hazard import general as haz_general
from openquake.engine.calculators.hazard.classical import (
    post_processing as cls_post_proc)
from openquake.engine.calculators.hazard.event_based import post_processing
from openquake.engine.db import models
from openquake.engine.input import logictree
from openquake.engine.utils import tasks


#: Always 1 for the computation of ground motion fields in the event-based
#: hazard calculator.
DEFAULT_GMF_REALIZATIONS = 1

# NB: beware of large caches; GmfRupture rows contain large arrays
inserter = writer.CacheInserter(models.GmfRupture, 10)


# Disabling pylint for 'Too many local variables'
# pylint: disable=R0914
@tasks.momotask
def compute_and_save_ses(task_mon, job_id, src_id, ses, seed):
    """
    Celery task for the stochastic event set calculator.

    Samples logic trees and calls the stochastic event set calculator.

    Once stochastic event sets are calculated, results will be saved to the
    database. See :class:`openquake.engine.db.models.SESCollection`.

    Optionally (specified in the job configuration using the
    `ground_motion_fields` parameter), GMFs can be computed from each rupture
    in each stochastic event set. GMFs are also saved to the database.

    :param int job_id:
        ID of the currently running job.
    :param src_id:
        Id of a parsed source model from which we will generate
        stochastic event sets/ruptures.
    :param lt_rlz:
        Logic Tree Realization object
    :param seed:
        Value for seeding numpy/scipy in the computation of stochastic event
        sets.
    """
    numpy.random.seed(seed)

    hc = models.HazardCalculation.objects.get(oqjob=job_id)
    lt_rlz = ses.ses_collection.lt_realization

    # complete_logic_tree_ses flag
    cmplt_lt_ses = None
    if hc.complete_logic_tree_ses:
        cmplt_lt_ses = models.SES.objects.get(
            ses_collection__output__oq_job=job_id,
            ordinal=None)

    # preparing sources
    ltp = logictree.LogicTreeProcessor(hc.id)
    apply_uncertainties = ltp.parse_source_model_logictree_path(
        lt_rlz.sm_lt_path)

    with task_mon.copy('reading source'):
        [src] = haz_general.gen_sources(
            [src_id], apply_uncertainties, hc.rupture_mesh_spacing,
            hc.width_of_mfd_bin, hc.area_source_discretization)

    src_filter = filters.source_site_distance_filter(hc.maximum_distance)
    rup_filter = filters.rupture_site_distance_filter(hc.maximum_distance)

    with task_mon.copy('reading site collection'):
        site_collection = hc.site_collection

    # Compute and save stochastic event sets
    # For each rupture generated, we can optionally calculate a GMF
    with task_mon.copy('computing ses ruptures'):
        # make copies of the hazardlib ruptures (which may contain
        # duplicates)
        ruptures = map(
            copy.copy,
            stochastic.stochastic_event_set_poissonian(
                [src], hc.investigation_time, site_collection,
                src_filter, rup_filter))
        # set the tag for each copy
        for i, r in enumerate(ruptures):
            r.tag = 'rlz=%02d|ses=%04d|src=%s|i=%03d' % (
                lt_rlz.ordinal, ses.ordinal, src.source_id, i)

    if ruptures:
        with task_mon.copy('saving ses ruptures'):
            _save_ses_ruptures(ses, ruptures, cmplt_lt_ses)


@tasks.momotask
def compute_and_save_gmf(task_mon, job_id, gsim, rupture):
    """
    """
    with task_mon.copy('reading calculation and site collection'):
        hc = models.HazardCalculation.objects.get(oqjob=job_id)
        imts = [haz_general.imt_to_hazardlib(x)
                for x in hc.intensity_measure_types]
        correl_model = None
        if hc.ground_motion_correlation_model is not None:
            correl_model = haz_general.get_correl_model(hc)
        site_collection = hc.site_collection

    with task_mon.copy('computing gmfs'):
        # Compute and save ground motion fields
        gmf_calc_kwargs = {
            'rupture': rupture.rupture,
            'sites': site_collection,
            'imts': imts,
            'gsim': gsim,
            'truncation_level': hc.truncation_level,
            'realizations': DEFAULT_GMF_REALIZATIONS,
            'correlation_model': correl_model,
            'rupture_site_filter': filters.rupture_site_distance_filter(
                hc.maximum_distance),
        }
        gmf_dict = gmf.ground_motion_fields(**gmf_calc_kwargs)

    with task_mon.copy('saving gmfs'):
        _save_gmf(rupture, gmf_dict, len(site_collection))


@transaction.commit_on_success(using='reslt_writer')
def _save_ses_ruptures(ses, ruptures, complete_logic_tree_ses):
    """
    Helper function for saving stochastic event set ruptures to the database.

    :param ses:
        A :class:`openquake.engine.db.models.SES` instance. This will be DB
        'container' for the new rupture record.
    :param rupture:
        A :class:`openquake.hazardlib.source.rupture.Rupture` instance.
    :param complete_logic_tree_ses:
        :class:`openquake.engine.db.models.SES` representing the `complete
        logic tree` stochastic event set.
        If not None, save a copy of the input `rupture` to this SES.
    """
    for r in ruptures:
        models.SESRupture.objects.create(ses=ses, rupture=r, tag=r.tag)

    if complete_logic_tree_ses is not None:
        for rupture in ruptures:
            models.SESRupture.objects.create(
                ses=complete_logic_tree_ses, rupture=rupture)


@transaction.commit_on_success(using='reslt_writer')
def _save_gmf(rupture, gmf_dict, nsites):
    """
    Helper method to save computed GMF data to the database.
    """
    for imt, gmvs in gmf_dict.iteritems():
        if all(gmv < 1E-6 for gmv in gmvs):
            return
        sa_period = None
        sa_damping = None
        if isinstance(imt, openquake.hazardlib.imt.SA):
            sa_period = imt.period
            sa_damping = imt.damping
        imt_name = imt.__class__.__name__
        inserter.add(
            models.GmfRupture(
                rupture=rupture,
                imt=imt_name,
                sa_period=sa_period,
                sa_damping=sa_damping,
                gmvs=list(gmvs.reshape(nsites))))
    inserter.flush()


class EventBasedHazardCalculator(haz_general.BaseHazardCalculator):
    """
    Probabilistic Event-Based hazard calculator. Computes stochastic event sets
    and (optionally) ground motion fields.
    """
    core_calc_task = compute_and_save_ses

    def task_arg_gen(self):
        """
        Loop through realizations and sources to generate a sequence of
        task arg tuples. Each tuple of args applies to a single task.
        Yielded results are tuples of the form job_id, src_ids, ses, task_seed
        (task_seed will be used to seed numpy for temporal occurence sampling).
        """
        hc = self.hc
        rnd = random.Random()
        rnd.seed(hc.random_seed)
        realizations = self._get_realizations()
        for lt_rlz in realizations:
            sources = models.SourceProgress.objects\
                .filter(is_complete=False, lt_realization=lt_rlz)\
                .order_by('id')\
                .values_list('parsed_source_id', flat=True)
            all_ses = models.SES.objects.filter(
                ses_collection__lt_realization=lt_rlz,
                ordinal__isnull=False).order_by('ordinal')
            for ses in all_ses:
                for src_id in sources:
                    seed = rnd.randint(0, models.MAX_SINT_32)
                    task_args = (self.job.id, src_id, ses, seed)
                    yield task_args

    def compute_and_save_gmf_arg_gen(self):
        """
        """
        for lt_rlz in self._get_realizations():
            ltp = logictree.LogicTreeProcessor(self.hc.id)
            gsims = ltp.parse_gmpe_logictree_path(lt_rlz.gsim_lt_path)
            ruptures = models.SESRupture.objects.filter(
                ses__ses_collection__lt_realization=lt_rlz)
            for rupture in ruptures:
                yield self.job.id, gsims[rupture.tectonic_region_type], rupture

    def execute(self):
        """
        Run compute_and_save_ses in parallel.
        """
        self.parallelize(self.core_calc_task, self.task_arg_gen())
        self.parallelize(compute_and_save_gmf,
                         self.compute_and_save_gmf_arg_gen())

    def initialize_ses_db_records(self, lt_rlz):
        """
        Create :class:`~openquake.engine.db.models.Output`,
        :class:`~openquake.engine.db.models.SESCollection` and
        :class:`~openquake.engine.db.models.SES` "container" records for
        a single realization.

        Stochastic event set ruptures computed for this realization will be
        associated to these containers.

        NOTE: Many tasks can contribute ruptures to the same SES.
        """
        output = models.Output.objects.create(
            owner=self.job.owner,
            oq_job=self.job,
            display_name='ses-coll-rlz-%s' % lt_rlz.id,
            output_type='ses')

        ses_coll = models.SESCollection.objects.create(
            output=output, lt_realization=lt_rlz)

        if self.job.hazard_calculation.ground_motion_fields:
            output = models.Output.objects.create(
                owner=self.job.owner,
                oq_job=self.job,
                display_name='gmf-rlz-%s' % lt_rlz.id,
                output_type='gmf')

            models.Gmf.objects.create(
                output=output, lt_realization=lt_rlz)

        all_ses = []
        for i in xrange(1, self.hc.ses_per_logic_tree_path + 1):
            all_ses.append(
                models.SES.objects.create(
                    ses_collection=ses_coll,
                    investigation_time=self.hc.investigation_time,
                    ordinal=i))
        return all_ses

    def initialize_complete_lt_ses_db_records(self):
        """
        Optional; if the user has requested to collect a `complete logic tree`
        stochastic event set (containing all ruptures from all realizations),
        initialize DB records for those results here.

        Throughout the course of the calculation, computed ruptures will be
        copied into this collection. See :func:`_save_ses_ruptures` for more
        info.
        """
        # `complete logic tree` SES
        clt_ses_output = models.Output.objects.create(
            owner=self.job.owner,
            oq_job=self.job,
            display_name='complete logic tree SES',
            output_type='complete_lt_ses')

        clt_ses_coll = models.SESCollection.objects.create(
            output=clt_ses_output)

        models.SES.objects.create(
            ses_collection=clt_ses_coll,
            investigation_time=self.hc.total_investigation_time())

        if self.hc.complete_logic_tree_gmf:
            clt_gmf_output = models.Output.objects.create(
                owner=self.job.owner,
                oq_job=self.job,
                display_name='complete logic tree GMF',
                output_type='complete_lt_gmf')
            models.Gmf.objects.create(output=clt_gmf_output)

    def pre_execute(self):
        """
        Do pre-execution work. At the moment, this work entails:
        parsing and initializing sources, parsing and initializing the
        site model (if there is one), parsing vulnerability and
        exposure files, and generating logic tree realizations. (The
        latter piece basically defines the work to be done in the
        `execute` phase.)
        """
        # Parse risk models.
        self.parse_risk_models()

        # Parse logic trees and create source Inputs.
        self.initialize_sources()

        # Deal with the site model and compute site data for the calculation
        # If no site model file was specified, reference parameters are used
        # for all sites.
        self.initialize_site_model()

        # Now bootstrap the logic tree realizations and related data.
        # This defines for us the "work" that needs to be done when we reach
        # the `execute` phase.
        rlz_callbacks = [self.initialize_ses_db_records]

        self.initialize_realizations(rlz_callbacks=rlz_callbacks)

        if self.job.hazard_calculation.complete_logic_tree_ses:
            self.initialize_complete_lt_ses_db_records()

        self.record_init_stats()

        num_sources = models.SourceProgress.objects.filter(
            is_complete=False,
            lt_realization__hazard_calculation=self.hc).count()
        self.progress['total'] = num_sources

        self.initialize_pr_data()

    def post_process(self):
        """
        If requested, perform additional processing of GMFs to produce hazard
        curves.
        """
        if self.hc.hazard_curves_from_gmfs:
            with self.monitor('generating hazard curves'):
                self.parallelize(
                    post_processing.gmf_to_hazard_curve_task,
                    post_processing.gmf_to_hazard_curve_arg_gen(self.job))

            # If `mean_hazard_curves` is True and/or `quantile_hazard_curves`
            # has some value (not an empty list), do this additional
            # post-processing.
            if self.hc.mean_hazard_curves or self.hc.quantile_hazard_curves:
                with self.monitor('generating mean/quantile curves'):
                    self.do_aggregate_post_proc()

            if self.hc.hazard_maps:
                with self.monitor('generating hazard maps'):
                    self.parallelize(
                        cls_post_proc.hazard_curves_to_hazard_map_task,
                        cls_post_proc.hazard_curves_to_hazard_map_task_arg_gen(
                            self.job))
