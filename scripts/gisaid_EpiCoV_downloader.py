#!/usr/bin/env python
"""

"""
__author__ = "Po-E Li, B10, LANL"
__copyright__ = "LANL 2020"
__license__ = "GPL"
__version__ = "1.0.0"
__email__ = "po-e@lanl.gov"

import os
import time
import sys
import argparse as ap
import json
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC


def parse_params():
    p = ap.ArgumentParser(prog='gisaid_EpiCoV_downloader.py',
                          description="""Download EpiCoV sequences from GISAID""")

    p.add_argument('-u', '--username',
                   metavar='[STR]', nargs=1, type=str, required=True,
                   help="GISAID username")

    p.add_argument('-p', '--password',
                   metavar='[STR]', nargs=1, type=str, required=True,
                   help="GISAID password")

    p.add_argument('-o', '--outdir',
                   metavar='[STR]', type=str, required=False, default=None,
                   help="Output directory")

    p.add_argument('-l', '--location',
                   metavar='[STR]', type=str, required=False, default=None,
                   help="sample location")

    p.add_argument('-cs', '--colstart',
                   metavar='[YYYY-MM-DD]', type=str, required=False, default=None,
                   help="collection starts date")

    p.add_argument('-ce', '--colend',
                   metavar='[YYYY-MM-DD]', type=str, required=False, default=None,
                   help="collection ends date")

    p.add_argument('-ss', '--substart',
                   metavar='[YYYY-MM-DD]', type=str, required=False, default=None,
                   help="submitssion starts date")

    p.add_argument('-se', '--subend',
                   metavar='[YYYY-MM-DD]', type=str, required=False, default=None,
                   help="submitssion ends date")

    p.add_argument('-cg', '--complete',
                   action='store_true', help='complete genome only')

    p.add_argument('-hc', '--highcoverage',
                   action='store_true', help='high coverage only')

    p.add_argument('-le', '--lowcoverageExcl',
                   action='store_true', help='low coverage excluding')

    p.add_argument('-t', '--timeout',
                   metavar='[INT]', type=int, required=False, default=90,
                   help="set action timeout seconds. Default is 90 secs.")

    p.add_argument('-r', '--retry',
                   metavar='[INT]', type=int, required=False, default=5,
                   help="retry how many times when the action fails. Default is 5 times.")

    p.add_argument('-i', '--interval',
                   metavar='[INT]', type=int, required=False, default=3,
                   help="time interval between retries in second(s). Default is 3 seconds.")

    p.add_argument('-m', '--meta',
                   action='store_true', help='download detail metadata (experimental, very slow)')

    p.add_argument('--normal',
                   action='store_true', help='run firefox in normal mode.')

    args_parsed = p.parse_args()
    if not args_parsed.outdir:
        args_parsed.outdir = os.getcwd()
    return args_parsed


def download_gisaid_EpiCoV(
        uname,     # username
        upass,     # password
        normal,  # normal mode (quite)
        wd,        # output dir
        loc,       # location
        cs,        # collection start date
        ce,        # collection end date
        ss,        # submission start date
        se,        # submission end date
        cg,        # complete genome only
        hc,        # high coverage only
        le,        # low coverage excluding
        to,        # timeout in sec
        rt,        # num of retry
        iv,        # interval in sec
        meta_dl    # also download meta
    ):
    """Download sequences and metadata from EpiCoV GISAID"""

    # output directory
    if not os.path.exists(wd):
        os.makedirs(wd, exist_ok=True)

    wd = os.path.abspath(wd)
    # GISAID_FASTA = f'{wd}/sequences.fasta.bz2'
    # GISAID_TABLE = f'{wd}/gisaid_cov2020_acknowledgement_table.xls'
    GISAID_DTL_JASON = f'{wd}/gisaid_detail_metadata.json'
    # GISAID_TSV   = f'{wd}/metadata.tsv.bz2'
    metadata = []

    # MIME types
    mime_types = "application/octet-stream"
    mime_types += ",application/excel,application/vnd.ms-excel"
    mime_types += ",application/pdf,application/x-pdf"
    mime_types += ",application/x-bzip2"
    mime_types += ",application/x-gzip,application/gzip"

    # start fresh
    try:
        os.remove(GISAID_DTL_JASON)
    except OSError:
        pass

    print("Opening browser...")
    profile = webdriver.FirefoxProfile()
    profile.set_preference("browser.download.folderList", 2)
    profile.set_preference("browser.download.manager.showWhenStarting", False)
    profile.set_preference("browser.download.dir", wd)
    profile.set_preference(
        "browser.helperApps.neverAsk.saveToDisk", mime_types)
    profile.set_preference(
        "plugin.disable_full_page_plugin_for_types", mime_types)
    profile.set_preference("pdfjs.disabled", True)

    options = Options()
    if not normal:
        options.add_argument("--headless")
    driver = webdriver.Firefox(firefox_profile=profile, options=options)

    # driverwait
    driver.implicitly_wait(20)
    wait = WebDriverWait(driver, to)

    # open GISAID
    print("Opening website GISAID...")
    driver.get('https://platform.gisaid.org/epi3/frontend')
    waiting_sys_timer(wait)
    print(driver.title)
    assert 'GISAID' in driver.title

    # login
    print("Logining to GISAID...")
    username = driver.find_element_by_name('login')
    username.send_keys(uname)
    password = driver.find_element_by_name('password')
    password.send_keys(upass)
    driver.execute_script("return doLogin();")

    waiting_sys_timer(wait)

    # navigate to EpiFlu
    print("Navigating to EpiCoV...")
    epicov_tab = driver.find_element_by_xpath("//div[@id='main_nav']//li[3]/a")
    epicov_tab.click()

    waiting_sys_timer(wait)

    # when user doesn't enter time/location, download all sequences and metadata
    if not (cs or ce or ss or se or loc):
        # download from downloads section
        print("Clicking downloads...")
        pd_button = wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//div[@class='sys-actionbar-bar']//div[3]")))
        pd_button.click()
        waiting_sys_timer(wait)

        # have to click the first row twice to start the iframe
        iframe = waiting_for_iframe(wait, driver, rt, iv)
        driver.switch_to.frame(iframe)
        waiting_sys_timer(wait)

        print("Downloading nextfasta...")
        dl_button = wait.until(EC.element_to_be_clickable(
            (By.XPATH, '//div[contains(text(), "nextfasta")]')))
        dl_button.click()
        waiting_sys_timer(wait)

        fn = wait_downloaded_filename(wait, driver, 180)
        print(f"Downloaded to {fn}.")

        waiting_sys_timer(wait)

        print("Downloading nextmeta...")
        dl_button = wait.until(EC.element_to_be_clickable(
            (By.XPATH, '//div[contains(text(), "nextmeta")]')))
        dl_button.click()

        fn = wait_downloaded_filename(wait, driver, 180)
        print(f"Downloaded to {fn}.")

        waiting_sys_timer(wait)

        # go back to main frame
        back_button = wait.until(EC.element_to_be_clickable(
            (By.XPATH, '//button[contains(text(), "Back")]')))
        back_button.click()

        driver.switch_to.default_content()
        waiting_sys_timer(wait)

    # have to reduce the range of genomes
    if cs or ce or ss or se or loc:
        print("Browsing EpiCoV...")
        browse_tab = wait.until(EC.element_to_be_clickable(
            (By.XPATH, '//*[contains(text(), "Browse")]')))
        browse_tab.click()
        waiting_sys_timer(wait)
        waiting_table_to_get_ready(wait)

        # set location
        if loc:
            print("Setting location...")
            loc_input = driver.find_element_by_xpath(
                "tr/td[contains(text(), 'Location')]/following-sibling::td/div/div/input"
            )
            loc_input.send_keys(loc)
            waiting_sys_timer(wait, 7)

        # set dates
        date_inputs = driver.find_elements_by_css_selector(
            "div.sys-form-fi-date input")
        dates = (cs, ce, ss, se)
        for dinput, date in zip(date_inputs, dates):
            if date:
                print("Setting date...")
                dinput.send_keys(date)

        ActionChains(driver).send_keys(Keys.ESCAPE).perform()
        waiting_sys_timer(wait, 7)

        # complete genome only
        if cg:
            print("complete genome only...")
            checkbox = driver.find_element_by_xpath('//input[@value="complete"]')
            checkbox.click()
            waiting_sys_timer(wait)

        # high coverage only
        if hc:
            print("high coverage only...")
            checkbox = driver.find_element_by_xpath('//input[@value="highq"]')
            checkbox.click()
            waiting_sys_timer(wait)

        # excluding low coverage
        if le:
            print("low coverage excluding...")
            checkbox = driver.find_element_by_xpath('//input[@value="lowco"]')
            checkbox.click()
            waiting_sys_timer(wait)

        # select all genomes
        print("Selecting all genomes...")
        button_sa = driver.find_element_by_css_selector("span.yui-dt-label input")
        button_sa.click()
        waiting_sys_timer(wait)

        # downloading sequence
        retry = 0
        while retry <= rt:
            try:
                print("Downloading sequences for selected genomes...")
                button = driver.find_element_by_xpath(
                    "//td[@class='sys-datatable-info']/button[contains(text(), 'Download')]")
                button.click()
                waiting_sys_timer(wait)

                # switch to iframe
                iframe = waiting_for_iframe(wait, driver, rt, iv)
                driver.switch_to.frame(iframe)
                waiting_sys_timer(wait)

                button = driver.find_element_by_xpath(
                    "//button[contains(text(), 'Download')]")
                button.click()
                waiting_sys_timer(wait)
                driver.switch_to.default_content()

                fn = wait_downloaded_filename(wait, driver, 180)
                print(f"Downloaded to {fn}.")

                break
            except:
                print(f"retrying...#{retry} in {iv} sec(s)")
                if retry == rt:
                    print("Unexpected error:", sys.exc_info())
                    sys.exit(1)
                else:
                    time.sleep(iv)
                    retry += 1

        # downloading metadata
        retry = 0
        while retry <= rt:
            try:
                print("Downloading acknowledgement table for selected genomes...")
                button = driver.find_element_by_xpath(
                    "//td[@class='sys-datatable-info']/button[contains(text(), 'Download')]")
                button.click()
                waiting_sys_timer(wait)

                # switch to iframe
                iframe = waiting_for_iframe(wait, driver, rt, iv)
                driver.switch_to.frame(iframe)
                waiting_sys_timer(wait)

                label = driver.find_element_by_xpath(
                    "//label[contains(text(), 'Acknowledgement Table')]")
                label.click()

                button = driver.find_element_by_xpath(
                    "//button[contains(text(), 'Download')]")
                button.click()

                waiting_sys_timer(wait)
                driver.switch_to.default_content()

                fn = wait_downloaded_filename(wait, driver, 180)
                print(f"Downloaded to {fn}.")

                break
            except:
                print(f"retrying...#{retry} in {iv} sec(s)")
                if retry == rt:
                    print("Unexpected error:", sys.exc_info())
                    sys.exit(1)
                else:
                    time.sleep(iv)
                    retry += 1

        # iterate each pages
        if meta_dl:
            page_num = 1
            print("Retrieving metadata...")
            while True:
                print(f"Starting processing page# {page_num}...")
                # retrieve tables
                tbody = wait.until(
                    EC.presence_of_element_located(
                        (By.XPATH, "//tbody[@class='yui-dt-data']"))
                )

                waiting_table_to_get_ready(wait)

                # interate each row
                for tr in tbody.find_elements_by_tag_name("tr"):
                    td = tr.find_element_by_tag_name("td")
                    driver.execute_script("arguments[0].scrollIntoView();", td)

                    # have to click the first row twice to start the iframe
                    iframe = None
                    record_elem = None
                    retry = 1
                    while retry <= rt:
                        try:
                            td.click()
                            waiting_sys_timer(wait)
                            iframe = driver.find_element_by_xpath("//iframe")
                            if iframe:
                                break
                            else:
                                raise
                        except:
                            print(f"retrying...#{retry} in {iv} sec(s)")
                            if retry == rt:
                                print("Failed")
                                sys.exit(1)
                            else:
                                time.sleep(iv)
                                retry += 1

                    driver.switch_to.frame(iframe)

                    # detect error: "An internal server error occurred."
                    # and "error-token: DYX47"
                    error_token = driver.find_element_by_xpath("//b")
                    if error_token:
                        error_token_text = error_token.text
                        if "error-token" in error_token.text:
                            print(
                                "[FATAL ERROR] A website internal server error occurred.")
                            print(error_token_text)
                            sys.exit(1)

                    # get the element of table with metadata
                    record_elem = wait.until(
                        EC.presence_of_element_located(
                            (By.XPATH, "//div[@class='packer']"))
                    )

                    # parse metadata
                    m = getMetadata(record_elem)
                    metadata.append(m)
                    print(f"{m['Accession ID']}\t{m['Virus name']}")

                    # get back
                    ActionChains(driver).send_keys(Keys.ESCAPE).perform()
                    time.sleep(1)
                    driver.switch_to.default_content()

                print(f"Compeleted page# {page_num}.")
                page_num += 1

                # go to the next page
                retry = 1
                button_next_page = None
                try:
                    button_next_page = driver.find_element_by_xpath(
                        f'//a[@page="{page_num}"]')
                except:
                    break

                if button_next_page:
                    print(f"Entering page# {page_num}...")
                    while retry <= rt:
                        try:
                            button_next_page.click()
                            time.sleep(10)
                            current_page = driver.find_element_by_xpath(
                                '//span[@class="yui-pg-current-page yui-pg-page"]').text
                            if current_page != str(page_num):
                                raise
                            else:
                                break
                        except:
                            print(f"retrying...#{retry} in {iv} sec(s)")
                            if retry == rt:
                                print("Failed")
                                sys.exit(1)
                            else:
                                time.sleep(iv)
                                retry += 1

            # writing metadata to JSON file
            print("Writing detail metadata...")
            with open(GISAID_DTL_JASON, 'w') as outfile:
                json.dump(metadata, outfile)

    # close driver
    driver.quit()


def getMetadata(record_elem):
    """parse out metadata from the table"""
    meta = {}
    table = record_elem.find_element_by_tag_name("table")
    last_attr = ""
    for tr in table.find_elements_by_tag_name("tr"):
        if tr.get_attribute("colspan") == "2":
            # skip titles
            continue
        else:
            tds = tr.find_elements_by_tag_name("td")
            if len(tds) == 2:
                attr = tds[0].text.strip(":")
                val = tds[1].text
                if attr == "Address":
                    attr = f"{last_attr} {attr.lower()}"
                    if attr == "Submission Date address":
                        attr = "Submitter address"
                meta[attr] = val
                last_attr = attr
    return meta


def waiting_sys_timer(wait, sec=1):
    """wait for system timer"""
    try:
        wait.until(EC.invisibility_of_element_located(
            (By.XPATH,  "//div[@id='sys_timer']")))
    except:
        pass
    time.sleep(sec)


def waiting_table_to_get_ready(wait, sec=1):
    """wait for the table to be loaded"""
    wait.until(EC.invisibility_of_element_located(
        (By.XPATH,  "//tbody[@class='yui-dt-message']")))
    time.sleep(sec)

def waiting_for_iframe(wait, driver, rt, iv):
    iframe = None
    retry = 1
    while retry <= rt:
        try:
            wait.until(EC.presence_of_element_located((By.XPATH, "//iframe")))
            iframe = driver.find_element_by_xpath("//iframe")
            if iframe:
                return iframe
            else:
                raise
        except:
            print(f"retrying...#{retry} in {iv} sec(s)")
            if retry == rt:
                print("Failed")
                sys.exit(1)
            else:
                time.sleep(iv)
                retry += 1

def wait_downloaded_filename(wait, driver, waitTime=180):
    driver.execute_script("window.open()")
    wait.until(EC.new_window_is_opened)
    driver.switch_to.window(driver.window_handles[-1])
    driver.get("about:downloads")
    time.sleep(1)

    endTime = time.time()+waitTime
    while True:
        try:
            progress = driver.execute_script("return document.querySelector('.downloadContainer progress:first-of-type').value")
            fileName = driver.execute_script("return document.querySelector('.downloadContainer description:first-of-type').value")

            while progress < 100:
                time.sleep(1)
                progress = driver.execute_script("return document.querySelector('.downloadContainer progress:first-of-type').value")

            driver.close()
            driver.switch_to.window(driver.window_handles[0])
            time.sleep(2)
            return fileName
        except:
            pass
        time.sleep(1)
        if time.time() > endTime:
            break


def main():

    argvs = parse_params()
    print(f"--- Ingest at {time.strftime('%Y-%m-%d %H:%M:%S')} ---")
    download_gisaid_EpiCoV(
        argvs.username,
        argvs.password,
        argvs.normal,
        argvs.outdir,
        argvs.location,
        argvs.colstart,
        argvs.colend,
        argvs.substart,
        argvs.subend,
        argvs.complete,
        argvs.highcoverage,
        argvs.lowcoverageExcl,
        argvs.timeout,
        argvs.retry,
        argvs.interval,
        argvs.meta
    )

if __name__ == "__main__":
    main()
