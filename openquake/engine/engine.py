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

"""Engine: A collection of fundamental functions for initializing and running
calculations."""

import ConfigParser
import csv
import getpass
import os
import sys
import warnings

import openquake.engine

from django.core import exceptions
from django import db as django_db
from lxml import etree

from openquake.engine import kvs
from openquake.engine import logs
from openquake.engine.db import models
from openquake.engine.job.validation import validate
from openquake.engine.supervising import supervisor
from openquake.engine.utils import monitor, get_calculator_class
from openquake.engine.writer import CacheInserter

from openquake import hazardlib
from openquake import risklib
from openquake import nrmllib


INPUT_TYPES = dict(models.INPUT_TYPE_CHOICES)
UNABLE_TO_DEL_HC_FMT = 'Unable to delete hazard calculation: %s'
UNABLE_TO_DEL_RC_FMT = 'Unable to delete risk calculation: %s'


def prepare_job(user_name="openquake", log_level='progress'):
    """
    Create job for the given user, return it.

    :param str username:
        Username of the user who owns/started this job. If the username doesn't
        exist, a user record for this name will be created.
    :param str log_level:
        Defaults to 'progress'. Specify a logging level for this job. This
        level can be passed, for example, from the command line interface using
        the `--log-level` directive.
    :returns:
        :class:`openquake.engine.db.models.OqJob` instance.
    """
    return models.OqJob.objects.create(
        user_name=user_name,
        log_level=log_level,
        oq_version=openquake.engine.__version__,
        nrml_version=nrmllib.__version__,
        hazardlib_version=hazardlib.__version__,
        risklib_version=risklib.__version__,
    )


def parse_config(source):
    """Parse a dictionary of parameters from an INI-style config file.

    :param source:
        File-like object containing the config parameters.
    :returns:
        Two dictionaries (as a 2-tuple). The first contains all of the
        parameters/values parsed from the job.ini file. The second contains
        absolute paths to all of the files referenced in the job.ini, keyed by
        the parameter name.
    """
    cp = ConfigParser.ConfigParser()
    cp.readfp(source)

    base_path = os.path.dirname(
        os.path.join(os.path.abspath('.'), source.name))
    params = dict(base_path=base_path)
    params['inputs'] = dict()

    # Directory containing the config file we're parsing.
    base_path = os.path.dirname(os.path.abspath(source.name))

    for sect in cp.sections():
        for key, value in cp.items(sect):
            if key == 'sites_csv':
                # Parse site coordinates from the csv file,
                # return as MULTIPOINT WKT:
                path = value
                if not os.path.isabs(path):
                    # It's a relative path
                    path = os.path.join(base_path, path)
                params['sites'] = _parse_sites_csv(open(path, 'r'))
            elif key.endswith('_file'):
                input_type = key[:-5]
                if not input_type in INPUT_TYPES:
                    raise ValueError(
                        'The parameter %s in the .ini file does '
                        'not correspond to a valid input type' % key)
                path = value
                # The `path` may be a path relative to the config file, or it
                # could be an absolute path.
                if not os.path.isabs(path):
                    # It's a relative path.
                    path = os.path.join(base_path, path)

                params['inputs'][input_type] = path
            else:
                params[key] = value

    # convert back to dict as defaultdict is not pickleable
    params['inputs'] = dict(params['inputs'])

    # load source inputs (the paths are the source_model_logic_tree)
    smlt = params['inputs'].get('source_model_logic_tree')
    if smlt:
        params['inputs']['source'] = [
            os.path.join(base_path, src_path)
            for src_path in _collect_source_model_paths(smlt)]

    return params


def _parse_sites_csv(fh):
    """
    Get sites of interest from a csv file. The csv file (``fh``) is expected to
    have 2 columns: lon,lat.

    :param fh:
        File-like containing lon,lat coordinates in csv format.

    :returns:
        MULTIPOINT WKT representing all of the sites in the csv file.
    """
    reader = csv.reader(fh)
    coords = []
    for lon, lat in reader:
        coords.append('%s %s' % (lon, lat))

    return 'MULTIPOINT(%s)' % ', '.join(coords)


def create_calculation(model, params):
    """
    Given a params `dict` parsed from the config file, create a
    :class:`~openquake.engine.db.models.HazardCalculation`.

    :param model:
        a Calculation class object
    :param dict params:
        Dictionary of parameter names and values. Parameter names should match
        exactly the field names of
        :class:`openquake.engine.db.model.HazardCalculation`.
    :returns:
        an instance of a newly created `model`
    """
    if "export_dir" in params:
        params["export_dir"] = os.path.abspath(params["export_dir"])

    calc_fields = model._meta.get_all_field_names()

    for param in set(params) - set(calc_fields):
        # FIXME(lp). Django 1.3 does not allow using _id fields in model
        # __init__. We will check these fields in pre-execute phase
        if param not in [
                'preloaded_exposure_model_id', 'hazard_output_id',
                'hazard_calculation_id']:
            msg = "Unknown parameter '%s'. Ignoring."
            msg %= param
            warnings.warn(msg, RuntimeWarning)
            params.pop(param)

    calc = model(**params)
    calc.full_clean()
    calc.save()

    return calc


def _collect_source_model_paths(smlt):
    """
    Given a path to a source model logic tree or a file-like, collect all of
    the soft-linked path names to the source models it contains and return them
    as a uniquified list (no duplicates).
    """
    src_paths = []
    tree = etree.parse(smlt)
    for branch_set in tree.xpath('//nrml:logicTreeBranchSet',
                                 namespaces=nrmllib.PARSE_NS_MAP):

        if branch_set.get('uncertaintyType') == 'sourceModel':
            for branch in branch_set.xpath(
                    './nrml:logicTreeBranch/nrml:uncertaintyModel',
                    namespaces=nrmllib.PARSE_NS_MAP):
                src_paths.append(branch.text)
    return sorted(list(set(src_paths)))


def _create_job_stats(job):
    """
    Helper function to create job stats, which implicitly records the start
    time for the job.

    :param job:
        :class:`openquake.engine.db.models.OqJob` instance.
    """
    models.JobStats.objects.create(oq_job=job)


# used by bin/openquake
def run_calc(job, log_level, log_file, exports, job_type, supervised=True):
    """
    Run a calculation.

    :param job:
        :class:`openquake.engine.db.model.OqJob` instance which references a
        valid :class:`openquake.engine.db.models.RiskCalculation` or
        :class:`openquake.engine.db.models.HazardCalculation`.
    :param str log_level:
        The desired logging level. Valid choices are 'debug', 'info',
        'progress', 'warn', 'error', and 'critical'.
    :param str log_file:
        Complete path (including file name) to file where logs will be written.
        If `None`, logging will just be printed to standard output.
    :param list exports:
        A (potentially empty) list of export targets. Currently only "xml" is
        supported.
    :param str job_type:
        'hazard' or 'risk'
    :param bool supervised:
        Defaults to `True`. If `True`, run OpenQuake with a supervisor process,
        which monitors the job executor process and collects log messages.
    """
    calc_mode = getattr(job, '%s_calculation' % job_type).calculation_mode
    calc = get_calculator_class(job_type, calc_mode)(job)

    # Create job stats, which implicitly records the start time for the job
    _create_job_stats(job)

    # Closing all db connections to make sure they're not shared between
    # supervisor and job executor processes.
    # Otherwise, if one of them closes the connection it immediately becomes
    # unavailable for others.
    if supervised:
        django_db.close_connection()

        job_pid = os.fork()

        if not job_pid:
            # calculation executor process
            try:
                _job_exec(job, log_level, exports, job_type, calc)
            except Exception, ex:
                logs.LOG.critical("Calculation failed with exception: '%s'"
                                  % str(ex))
                raise
            finally:
                job.is_running = False
                job.save()
            return

        supervisor_pid = os.fork()
        if not supervisor_pid:
            # supervisor process
            logs.set_logger_level(logs.logging.root, log_level)
            # TODO: deal with KVS garbage collection
            supervisor.supervise(job_pid, job.id, log_file=log_file)
            return

        # parent process

        # ignore Ctrl-C as well as supervisor process does. thus only
        # job executor terminates on SIGINT
        supervisor.ignore_sigint()
        # wait till both child processes are done
        os.waitpid(job_pid, 0)
        os.waitpid(supervisor_pid, 0)
    else:
        try:
            _job_exec(job, log_level, exports, job_type, calc)
        except Exception, ex:
            logs.LOG.critical("Calculation failed with exception: '%s'"
                              % str(ex))
            raise
        finally:
            job.is_running = False
            job.save()
            # Normally the supervisor process does this, but since we don't
            # have one in this case, we have to call the cleanup manually.
            supervisor.cleanup_after_job(job.id)

    # Refresh the job record, in case we are forking and another process has
    # modified the job state.
    return _get_job(job.id)


def _get_job(job_id):
    """
    Helper function to get a job object by ID. Makes testing/mocking easier.
    """
    return models.OqJob.objects.get(id=job_id)


def _job_exec(job, log_level, exports, job_type, calc):
    """
    Abstraction of some general job execution procedures.

    Parameters are the same as :func:`run_calc`, except for ``supervised``
    which is not included. Also ``calc`` is an instance of the calculator class
    which is passed to :func:`_do_run_calc`.
    """
    logs.init_logs_amqp_send(level=log_level, calc_domain=job_type,
                             calc_id=job.calculation.id)
    # run the job
    job.is_running = True
    job.save()
    kvs.mark_job_as_current(job.id)
    _do_run_calc(job, exports, calc, job_type)


def _switch_to_job_phase(job, ctype, status):
    """Switch to a particular phase of execution.

    This involves creating a `job_phase_stats` record and logging the new
    status.

    :param job:
        An :class:`~openquake.engine.db.models.OqJob` instance.
    :param str ctype: calculation type (hazard|risk)
    :param str status: one of the following: pre_executing, executing,
        post_executing, post_processing, export, clean_up, complete
    """
    job.status = status
    job.save()
    models.JobPhaseStats.objects.create(oq_job=job, job_status=status,
                                        ctype=ctype)
    logs.LOG.progress("%s (%s)" % (status, ctype))
    if status == "executing" and not openquake.engine.no_distribute():
        # Record the compute nodes that were available at the beginning of the
        # execute phase so we can detect failed nodes later.
        failed_nodes = monitor.count_failed_nodes(job)
        if failed_nodes == -1:
            logs.LOG.critical("No live compute nodes, aborting calculation")
            sys.exit(1)


def _do_run_calc(job, exports, calc, job_type):
    """
    Step through all of the phases of a calculation, updating the job
    status at each phase.

    :param job:
        An :class:`~openquake.engine.db.models.OqJob` instance.
    :param list exports:
        a (potentially empty) list of export targets, currently only "xml" is
        supported
    :returns:
        The input job object when the calculation completes.
    """
    # - Run the calculation
    _switch_to_job_phase(job, job_type, "pre_executing")

    calc.pre_execute()

    _switch_to_job_phase(job, job_type, "executing")
    calc.execute()

    _switch_to_job_phase(job, job_type, "post_executing")
    calc.post_execute()

    _switch_to_job_phase(job, job_type, "post_processing")
    calc.post_process()

    _switch_to_job_phase(job, job_type, "export")
    calc.export(exports=exports)

    _switch_to_job_phase(job, job_type, "clean_up")
    calc.clean_up()

    CacheInserter.flushall()  # flush caches into the db

    _switch_to_job_phase(job, job_type, "complete")
    logs.LOG.debug("*> complete")

    return job


def del_haz_calc(hc_id):
    """
    Delete a hazard calculation and all associated outputs.

    :param hc_id:
        ID of a :class:`~openquake.engine.db.models.HazardCalculation`.
    """
    try:
        hc = models.HazardCalculation.objects.get(id=hc_id)
    except exceptions.ObjectDoesNotExist:
        raise RuntimeError('Unable to delete hazard calculation: '
                           'ID=%s does not exist' % hc_id)

    user = getpass.getuser()
    if hc.oqjob.user_name == user:
        # we are allowed to delete this

        # but first, check if any risk calculations are referencing any of our
        # outputs, or the hazard calculation itself
        msg = UNABLE_TO_DEL_HC_FMT % (
            'The following risk calculations are referencing this hazard'
            ' calculation: %s'
        )

        # check for a reference to hazard outputs
        assoc_outputs = models.RiskCalculation.objects.filter(
            hazard_output__oq_job__hazard_calculation=hc_id
        )
        if assoc_outputs.count() > 0:
            raise RuntimeError(msg % ', '.join([str(x.id)
                                                for x in assoc_outputs]))

        # check for a reference to the hazard calculation itself
        assoc_calcs = models.RiskCalculation.objects.filter(
            hazard_calculation=hc_id
        )
        if assoc_calcs.count() > 0:
            raise RuntimeError(msg % ', '.join([str(x.id)
                                                for x in assoc_calcs]))

        # No risk calculation are referencing what we want to delete.
        # Carry on with the deletion.
        hc.delete(using='admin')
    else:
        # this doesn't belong to the current user
        raise RuntimeError(UNABLE_TO_DEL_HC_FMT % 'Access denied')


def del_risk_calc(rc_id):
    """
    Delete a risk calculation and all associated outputs.

    :param hc_id:
        ID of a :class:`~openquake.engine.db.models.RiskCalculation`.
    """
    try:
        rc = models.RiskCalculation.objects.get(id=rc_id)
    except exceptions.ObjectDoesNotExist:
        raise RuntimeError('Unable to delete risk calculation: '
                           'ID=%s does not exist' % rc_id)

    if rc.oqjob.user_name == getpass.getuser():
        # we are allowed to delete this
        rc.delete(using='admin')
    else:
        # this doesn't belong to the current user
        raise RuntimeError('Unable to delete risk calculation: '
                           'Access denied')


def run_hazard(cfg_file, log_level, log_file, exports):
    """
    Run a hazard calculation using the specified config file and other options.

    :param str cfg_file:
        Path to calculation config (INI-style) file.
    :param str log_level:
        'debug', 'info', 'warn', 'error', or 'critical'
    :param str log_file:
        Path to log file.
    :param list exports:
        A list of export types requested by the user. Currently only 'xml'
        is supported.
    """
    try:
        if log_file is not None:
            touch_log_file(log_file)

        job = haz_job_from_file(
            cfg_file, getpass.getuser(), log_level, exports
        )

        # Initialize the supervisor, instantiate the calculator,
        # and run the calculation.
        completed_job = run_calc(
            job, log_level, log_file, exports, 'hazard'
        )
        if completed_job is not None:
            # We check for `None` here because the supervisor and executor
            # process forks return to here as well. We want to ignore them.
            if completed_job.status == 'complete':
                print 'Calculation %d results:' % (
                    completed_job.hazard_calculation.id)
                list_hazard_outputs(completed_job.hazard_calculation.id)
            else:
                complain_and_exit('Calculation %d failed'
                                  % completed_job.hazard_calculation.id,
                                  exit_code=1)
    except IOError as e:
        print str(e)
    except Exception as e:
        raise


@django_db.transaction.commit_on_success
def haz_job_from_file(cfg_file_path, username, log_level, exports):
    """
    Create a full hazard job profile from a job config file.

    :param str cfg_file_path:
        Path to the job.ini.
    :param str username:
        The user who will own this job profile and all results.
    :param str log_level:
        Desired log level.
    :param exports:
        List of desired export types.

    :returns:
        :class:`openquake.engine.db.models.OqJob` object
    :raises:
        `RuntimeError` if the input job configuration is not valid
    """
    # create the job
    job = prepare_job(user_name=username, log_level=log_level)

    # read calculation params and create the calculation profile
    params = parse_config(open(cfg_file_path, 'r'))
    job.hazard_calculation = create_calculation(
        models.HazardCalculation, params)
    job.save()

    # validate and raise if there are any problems
    error_message = validate(job, 'hazard', params, exports)
    if error_message:
        raise RuntimeError(error_message)

    return job


def list_hazard_outputs(hc_id):
    """
    List the outputs for a given
    :class:`~openquake.engine.db.models.HazardCalculation`.

    :param hc_id:
        ID of a hazard calculation.
    """
    print_outputs_summary(get_hazard_outputs(hc_id))


def get_hazard_outputs(hc_id):
    """
    :param hc_id:
        ID of a hazard calculation.
    :returns:
        A sequence of :class:`openquake.engine.db.models.Output` objects
    """
    return models.Output.objects.filter(oq_job__hazard_calculation=hc_id)


def touch_log_file(log_file):
    """
    If a log file destination is specified, attempt to open the file in
    'append' mode ('a'). If the specified file is not writable, an
    :exc:`IOError` will be raised.
    """
    try:
        open(os.path.abspath(log_file), 'a').close()
    except IOError as e:
        raise IOError('Error writing to log file %s: %s'
                      % (log_file, e.strerror))


def complain_and_exit(msg, exit_code=0):
    """
    Print a ``msg`` and exit the current process with the given ``exit_code``.
    """
    print msg
    sys.exit(exit_code)


def print_outputs_summary(outputs):
    """
    List of :class:`openquake.engine.db.models.Output` objects.
    """
    if len(outputs) > 0:
        print 'id | output_type | name'
        for o in outputs.order_by('output_type'):
            print '%s | %s | %s' % (
                o.id, o.get_output_type_display(), o.display_name)


def run_risk(cfg_file, log_level, log_file, exports, hazard_output_id=None,
             hazard_calculation_id=None):
    """
    Run a risk calculation using the specified config file and other options.
    One of hazard_output_id or hazard_calculation_id must be specified.

    :param str cfg_file:
        Path to calculation config (INI-style) file.
    :param str log_level:
        'debug', 'info', 'warn', 'error', or 'critical'
    :param str log_file:
        Path to log file.
    :param list exports:
        A list of export types requested by the user. Currently only 'xml'
        is supported.
    :param str hazard_ouput_id:
        The Hazard Output ID used by the risk calculation (can be None)
    :param str hazard_calculation_id:
        The Hazard Calculation ID used by the risk calculation (can be None)
    """
    try:
        if log_file is not None:
            touch_log_file(log_file)

        job = risk_job_from_file(
            cfg_file, getpass.getuser(), log_level, exports, hazard_output_id,
            hazard_calculation_id
        )

        # Initialize the supervisor, instantiate the calculator,
        # and run the calculation.
        completed_job = run_calc(
            job, log_level, log_file, exports, 'risk'
        )
        if completed_job is not None:
            # We check for `None` here because the supervisor and executor
            # process forks return to here as well. We want to ignore them.
            if completed_job.status == 'complete':
                print 'Calculation %d results:' % (
                    completed_job.risk_calculation.id)
                list_risk_outputs(completed_job.risk_calculation.id)
            else:
                complain_and_exit('Calculation %s failed'
                                  % completed_job.risk_calculation.id,
                                  exit_code=1)
    except IOError as e:
        print str(e)
    except Exception as e:
        raise


@django_db.transaction.commit_on_success
def risk_job_from_file(cfg_file_path, username, log_level, exports,
                       hazard_output_id=None, hazard_calculation_id=None):
    """
    Create a full risk job profile from a job config file.

    :param str cfg_file_path:
        Path to the job.ini.
    :param str username:
        The user who will own this job profile and all results.
    :param str log_level:
        Desired log level.
    :param exports:
        List of desired export types.
    :param int hazard_output_id:
        ID of a hazard output to use as input to this calculation. Specify
        this xor ``hazard_calculation_id``.
    :param int hazard_calculation_id:
        ID of a complete hazard calculation to use as input to this
        calculation. Specify this xor ``hazard_output_id``.

    :returns:
        :class:`openquake.engine.db.models.OqJob` object
    :raises:
        `RuntimeError` if the input job configuration is not valid
    """
    assert not(hazard_output_id is None and hazard_calculation_id is None), (
        "Must specify either `hazard_output_id` or `hazard_calculation_id`, "
        "and not both"
    )
    # create the job
    job = prepare_job(user_name=username, log_level=log_level)

    # read calculation params and create the calculation profile
    params = parse_config(open(cfg_file_path, 'r'))
    # Add the hazard output id to the risk calculation constructor args
    params.update(dict(hazard_output_id=hazard_output_id,
                       hazard_calculation_id=hazard_calculation_id))

    calculation = create_calculation(models.RiskCalculation, params)
    job.risk_calculation = calculation
    job.save()

    error_message = validate(job, 'risk', params,  exports)
    if error_message:
        raise RuntimeError(error_message)

    return job


def list_risk_outputs(rc_id):
    """
    List the outputs for a given
    :class:`~openquake.engine.db.models.RiskCalculation`.

    :param rc_id:
        ID of a risk calculation.
    """
    print_outputs_summary(get_risk_outputs(rc_id))


def get_risk_outputs(rc_id):
    """
    :param rc_id:
        ID of a risk calculation.
    :returns:
        A sequence of :class:`openquake.engine.db.models.Output` objects
    """
    return models.Output.objects.filter(oq_job__risk_calculation=rc_id)
