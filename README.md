
.. -*- mode: rst -*-

.. <img align="right" width="100" height="100" src="https://github.com/evantwist/neuromonitoring_dashboard/assets/169645691/92209640-884b-4c41-8666-d989584ea7c4">

**PANDA** is an open-source *Pediatric Algorithm for Neuromonitoring and Dynamic Autoregulation* written in MATLAB (version 9.13) using the Statistics and Machine Learning Toolbox. Developed at the Erasmus MC Sophia Children's Hospital pediatric intensive care unit (PICU) in collaboration with Delft University of Technology. May 2024.

The main functions of PANDA are:

* Calculation of the pressure-reactivity index (PRx) from ICP and MAP data
* Calculation of the optimal CPP (CPPopt) and optimal range
* Calculation of ICP insults and color-coded visualization of correlation with outcome

PANDA is designed for users who want to perform neuromonitoring in pediatric traumatic brain injury with regard to cerebral autoregulation and obtain optimal targets for PRx, CPP and ICP that were shown to be associated with outcome at one year post-injury. The target values can be integrated into a comprehensive neuromonitoring dashboard.

**What are the prerequisites for using PANDA?**

To use PANDA all you need is:

- Some basic knowledge of MATLAB, especially working with arrays
- Some ICP and MAP data and optionally CPP data; ideally formatted as synchronous time series.

**I have ICP and MAP time series, how can I get started with using PANDA?**

If you have synchronized time series of ICP and MAP you can use the Import Tool and/or automatically generate scripts from here, please refer to `Matlab Documentation <https://nl.mathworks.com/help/matlab/import_export/import-data-interactively.html>`_ or for .csv/.xlsx/.txt/.mat-files you can use readmatrix(filename.csv)/readtable(filename.xlsx)/load(filename.mat)/fscanf(filename).

The next step is to clean the raw data of artefacts and applying filtering to remove noise from pulse and respiration. A simple preprocessing pipeline is shown below:


.. code-block:: MATLAB

  proc_signal = struct();
  % Read the file
  raw = readtable('raw_ICP_signal.xlsx');
  % Define physiological range for artefact detection
  sudden_deflections = 0.25;        % Percentage in or decrease from adjacent samples
  range_ll = 30;                    % mmHg
  range_ul = 160;                   % mmHg
  % Apply in loop
  for iii = 11:length(raw)
    if raw(ii) <= raw_signal(ii-1)*(1+sudden_deflections) || raw(ii) <= raw_signal(ii-1)*(1-sudden_deflections) || ...
    raw(ii) > range_ul || raw(ii) < range_ll
      proc_signal.clean(ii) = NaN;
      proc_signal.clean(ii-10:ii+60) = NaN; 
    elseif length(proc_signal) > ii
      proc_signal(ii) = NaN;
    else
      proc_signal(ii) = raw(ii);
    end
  end

  proc_signal = fillmissing(proc_signal, 'movmean', 100);
    

Afterwards you can dive in and attempt to:

* 'PRx calculation': calculate PRx  from your filtered ICP and MAP data (func_PRx).
* 'CPPopt calculation': use the acquired PRx alongside your ICP data to derive CPPopt target and optimal range (func_CPPopt).
* 'ICP insult analysis': analyze ICP on an insult basis, based on its intensity and duration, and correlate the average number of ICP insults with outcome in a color coded plot (func_ICP).


Development
~~~~~~~~~~~
PANDA was created in collaboration between the PICU of Erasmus MC Sophia Children's Hospital, Rotterdam, The Netherlands in collaboration with the TU Delft, Delft, The Netherlands. The initial project was created by Naomi Kathanarathan, Jan Willem Kuiper and Rogier de Jonge, pediatric-intensivists at the `Erasmus MC PICU <https://www.erasmusmc.nl/nl-nl/sophia/patientenzorg/specialismen/intensive-care-kinderen>`_ and `Alfred Schouten <https://www.tudelft.nl/staff/a.c.schouten/>`_, professor (Bio-)mechnical engineering at the TU Delft. The initial syntax was written by Bart Formsma and Tahisa Robles under supervision of Eris van Twist, and is maintained by `Eris van Twist <https://www.linkedin.com/in/eris-van-twist-423609183/>`_, a clinical technologist and researcher at Erasmus MC.

Note that this algorithm is provided with NO WARRANTY OF ANY KIND. Findings have been validated using retrospective data of the Erasmus MC PICU and that the algorithm has not been validated for real-time clinical use.

Feel free to contact to make a contribution, open an issue or submit a pull request.

Citation
~~~~~~~~

To cite PANDA, please refer to this GitHub as article approval is still pending.
