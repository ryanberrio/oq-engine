/*
    Copyright (c) 2010-2011, GEM Foundation.

    OpenQuake is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License version 3
    only, as published by the Free Software Foundation.

    OpenQuake is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License version 3 for more details
    (a copy is included in the LICENSE file that accompanied this code).

    You should have received a copy of the GNU Lesser General Public License
    version 3 along with OpenQuake.  If not, see
    <http://www.gnu.org/licenses/lgpl-3.0.txt> for a copy of the LGPLv3 License.
*/

-- Fix the status value constraints for uiapi.upload and uiapi.oq_job

ALTER TABLE uiapi.upload DROP CONSTRAINT upload_status_value;
ALTER TABLE uiapi.upload ADD CONSTRAINT upload_status_value
        CHECK(status IN ('pending', 'running', 'failed', 'succeeded'));

ALTER TABLE uiapi.oq_job DROP CONSTRAINT job_status_value;
ALTER TABLE uiapi.oq_job ADD CONSTRAINT job_status_value
        CHECK(status IN ('pending', 'running', 'failed', 'succeeded'));