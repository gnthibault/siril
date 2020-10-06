# Contributing to Siril

First, thanks for taking the time to contribute!

## How Can I Contribute?

Siril is Free Software and a part of the GNU Project and you are welcome to contribute to this project. There are many ways to do it.

 * program new features,
 * report bugs (errors in the program),
 * test existing features and provide feedback,
 * add and improve documentation,
 * translate Siril to your own language,
 * translate the documentation,
 * donate
 
### Getting last version

In order to try the last development version we recommend to compile the sources (see README). Nevertheless it is possible to get nightly builds for Linux and Windows by following these direct links:

 * GNU/Linux (x86_64): <https://gitlab.com/free-astro/siril/-/jobs/artifacts/master/download?job=appimage-nightly>
 * Windows (64bits): <https://gitlab.com/free-astro/siril/-/jobs/artifacts/master/download?job=win64-nightly>
 
 **Test builds are for testing purpose only. They have not been human-tested, it relies on regularly modified development code. So please do not use it for production!**
 
### Reporting Bugs

Reporting the bugs that you will encounter is very important to the development, it helps the developers to make Siril more stable and more bug free. If you have some programming skills you can attach a patch to your bug report, we will be happy to apply it.

#### Before Submitting A Bug Report

Before creating bug reports, please check [this list](https://gitlab.com/free-astro/siril/issues) as you might find out that you don't need to create one. Also, make sure you are using last stable version or last git version before submitting the ticket. 
When you are creating a bug report, please include as many details as possible. Fill out [the required template](https://gitlab.com/free-astro/siril/blob/master/.gitlab/issue_templates/bug.md), the information it asks for helps us resolve issues faster.

### Translation

We are looking for volunteer translators, for the software and for the documentation. No programming experience is required; you just need to download [ Poedit ](https://poedit.net/) for your OS, and generate *.pot and *.po files by making:

    ninja siril-pot -C _build
    ninja siril-update-po -C _build
    
Once done, you can either open the *.pot file for a new translation or the *.po file of the language you want to contribute to an existing one. Then, you can open a merge request or a bug report and attach the new *.po file. We do not need *.mo that are compiled files.

It is also possible to help translating the documentation. To do that, feel free to open a new ticket.

### Donate

If you like the software, please help us by contributing with the [ Donate button ](https://www.siril.org/#support-us) of the website. Siril takes us a lot of time and we still have to pay for the servers. 
